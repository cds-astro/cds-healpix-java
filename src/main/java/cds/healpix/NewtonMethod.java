// Copyright 2017-2018 - Universite de Strasbourg/CNRS
// The CDS HEALPix library is developped by the Centre de Donnees
// astronomiques de Strasbourgs (CDS) from the following external papers:
//  - [Gorsky2005]     - "HEALPix: A Framework for High-Resolution Discretization and
//                       Fast Analysis of Data Distributed on the Sphere"
//                       http://adsabs.harvard.edu/abs/2005ApJ...622..759G
//  - [Calabretta2004] - "Mapping on the HEALPix grid"
//                       http://adsabs.harvard.edu/abs/2004astro.ph.12607C
//  - [Calabretta2007] - "Mapping on the HEALPix grid"
//                       http://adsabs.harvard.edu/abs/2007MNRAS.381..865C
//  - [Reinecke2015]   - "Efficient data structures for masks on 2D grids"
//                       http://adsabs.harvard.edu/abs/2015A&A...580A.132R
// It is distributed under the terms of the BSD License 2.0
//
// This file is part of the CDS HEALPix library.
//

package cds.healpix;


import java.io.File;
import java.io.IOException;
// import java.nio.file.Files;
// import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import cds.healpix.common.sphgeom.CooXYZ;
import cds.healpix.common.sphgeom.Vect3D;

// import org.slf4j.Logger;
// import org.slf4j.LoggerFactory;

import static java.lang.Math.PI;

import static cds.healpix.common.sphgeom.CooXYZ.crossProd;

import static cds.healpix.Healpix.TRANSITION_Z;
import static cds.healpix.common.math.Math.FOUR_OVER_PI;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.sqrt;
import static cds.healpix.common.math.Math.abs;
import static cds.healpix.common.math.Math.tan;

final class NewtonMethod {

  //  private static final Logger LOG = LoggerFactory.getLogger("fr.unistra.cds.healpix.cone");

  private static final double SQRT_2_OVER_3 = 0.81649658092772603272;

  
  
  public static List<CooXYZ> arcSpecialPoints(CooXYZ p1, CooXYZ p2, double zEpsMax, final int nIterMax) {
     // Ensure p1.z() < p2.z()
    if (p1.z() > p2.z()) {
      final CooXYZ tmp = p1;
      p1 = p2;
      p2 = tmp;
    }
    final List<CooXYZ> res = new ArrayList<CooXYZ>();
    if (TRANSITION_Z <= p1.z() || p2.z() <= -TRANSITION_Z) { // NPC only or SPC only
      final CooXYZ r = arcSpecialPointInPc(p1, p2, zEpsMax, nIterMax);
      if (r != null) {
        res.add(r);
      } 
    } else if (-TRANSITION_Z <= p1.z() && p2.z() <= TRANSITION_Z) { // EQR only
      final CooXYZ r = arcSpecialPointInEqr(p1, p2, zEpsMax, nIterMax);
      if (r != null) {
        res.add(r);
      }
    } else if (p1.z() < -TRANSITION_Z) { // SPC, EQR (and maybe NPC)
      final CooXYZ vEqrSouth = new CooXYZ(intersectWithTransitionLatSpc(p1, p2));
      CooXYZ r = arcSpecialPointInPc(p1, vEqrSouth, zEpsMax, nIterMax);
      if (r != null) {
        res.add(r);
      }
      if (p2.z() <= TRANSITION_Z) {
       r = arcSpecialPointInEqr(vEqrSouth, p2, zEpsMax, nIterMax);
        if (r != null) {
          res.add(r);
        }
      } else {
        final CooXYZ vEqrNorth = new CooXYZ(intersectWithTransitionLatNpc(p1, p2));
        r = arcSpecialPointInEqr(vEqrSouth, vEqrNorth, zEpsMax, nIterMax);
        if (r != null) {
          res.add(r);
        }
        r = arcSpecialPointInPc(vEqrNorth, p2, zEpsMax, nIterMax);
        if (r != null) {
          res.add(r);
        }
      }
    } else { // both EQR and NPC
      final CooXYZ vEqrNorth = new CooXYZ(intersectWithTransitionLatNpc(p1, p2));
      CooXYZ r = arcSpecialPointInEqr(p1, vEqrNorth, zEpsMax, nIterMax);
      if (r != null) {
        res.add(r);
      }
      r = arcSpecialPointInPc(vEqrNorth, p2, zEpsMax, nIterMax);
      if (r != null) {
        res.add(r);
      }
    }
    return res;
  }
  
  /**
   * Returns the coordinate Y (=3/2 sin(lat)) of the point on the cone such that d(DeltaX)/dY = +-1,
   * i.e. the point on the cone having the equation Y = (+-)X + b with the larger (smallest)
   * possible y-intercept (b).
   * The function uses the Newton-Raphson method to solve (DeltaX)/dY -+ 1 = 0.
   * @param yStart the starting point of the Newton-Raphson method. A goog choise is
   *        3/2 * sin(lat + rRad * 0.8) for positive slopes, and
   *        3/2 * sin(lat - rRad * 0.8) for negative slopes
   * @param sinOfConeCenterLat the sine of the latitude of the center of the cone (= z_0)
   * @param coneRadiusRad the raidus of the cone (in radians)
   * @param positiveSlope {@code true} if we look for lines of equations Y = X + b,
   *                      {@code false} if we look for lines of equations Y = -X + b,
   *                       
   * @return
   */
  /*public static double newtonSolveEquatorialZone(double yStart,
      double sinOfConeCenterLat, double coneRadiusRad, boolean positiveSlope,
      double yEpsMax, int nIterMax) { // yEpsMax = 10e-8 => angle of ~= 2 mas
    final double a = positiveSlope ? FOUR_OVER_PI : - FOUR_OVER_PI;    
    final double sqrtW0 = Math.sqrt(1 - sinOfConeCenterLat * sinOfConeCenterLat);
    final double u0 = (1 - 0.5 * coneRadiusRad * coneRadiusRad) / sqrtW0;
    final double v0 = Math.abs(sinOfConeCenterLat) / sqrtW0;

    double prevY = -1, y = yStart;

    for (int i = 0; i < nIterMax && Math.abs(y - prevY) > yEpsMax; i++) {
      final double sqrt1mW =  Math.abs(y) * TRANSITION_Z;
      final double w = 1 - sqrt1mW * sqrt1mW;
      final double p = v0 - u0 * sqrt1mW;
      final double m = u0 - v0 * sqrt1mW;
      final double o = Math.sqrt(w - m * m);z0 - q * n
      final double f = Math.copySign(p / (w * o), y) - a;
      final double df = (2 * sqrt1mW * p) / (w * w * o) - u0 / (w * o)
        + (p / (w * o * o * o)) * (sqrt1mW - v0 * m);
      prevY = y;
      y += f / df;
      System.out.println("y: " + Math.toDegrees(Math.asin(prevY / 1.5)) + "; f: " + f + "; df: " + df);
    }
    return y;
  }*/

  /**
   * Returns the coordinate z (=sin(lat)) of the point on the cone such that d(DeltaX)/dY = +-1,
   * i.e. the point on the cone having the equation Y = (+-)X + b with the larger (smallest)
   * possible y-intercept (b).
   * 
   * Ok for radius > 20 mas, then if radius < 20 mas = 1e-7 rad
   * then
   * 1 - 2*sin^2(r/2) = 1 - 1e-14 which is close from the precision of a double
   * 
   * The function uses the Newton-Raphson method to solve (DeltaX)/dY -+ 1 = 0.
   * @param zStart the starting point of the Newton-Raphson method. A good choice is
   *        3/2 * sin(lat - 0.5 * rRad) for positive slopes, and
   *        3/2 * sin(lat + 0.5 * rRad) for negative slopes
   * @param sinOfConeCenterLat the sine of the latitude of the center of the cone (= z_0)
   * @param coneRadiusRad the radius of the cone (in radians)
   * @param positiveSlope {@code true} if we look for lines of equations Y = X + b 
   *                      (i.e. south part of the cone),
   *                      {@code false} if we look for lines of equations Y = -X + b
   *                      (e.e. north part of the cone)
   * @param zEpsMax the precision on z to stop the iteration (choose e.g. coneRadiusRad/1000).
   * @param nIterMax the maximum number of iteration in the Newton-Raphson algorithm
   * @return
   */
  public static double newtonSolveEquatorialZone(double zStart,
      double sinOfConeCenterLat, double twoSineOfHalfConeRadius, boolean positiveSlope,
      double zEpsMax, int nIterMax) { // yEpsMax = 10e-8 => angle of ~= 2 mas
    final double w0 = 1 - sinOfConeCenterLat * sinOfConeCenterLat;
    final double r = (1 - 0.5 * twoSineOfHalfConeRadius * twoSineOfHalfConeRadius);
    final double z0 = sinOfConeCenterLat;
    final double cte = 1.5 * (positiveSlope ? PI_OVER_FOUR : -PI_OVER_FOUR); // 1.5 = 1 / TRANSITION_Z

    double zEps = 1, z = zStart;
    for (int i = 0; i < nIterMax && abs(zEps) > zEpsMax; i++) {
      zEps = fOverDfEqr(z, z0, w0, cte, r);
      z -= zEps;
    }
    return z;
  }

  /**
   * Check if the great-circle arc (defined by the smallest distance between the two given points)
   * contains a 'special point', i.e. a point such such that a tangent line arc on the projection 
   * plane has a slope equals to `+-1`, i.e. `d(DeltaX)/dY = +-1`.
   * @param p1 first point of the great-circle arc
   * @param p2 second point of the great_-circle arc
   * @param zEpsMax the target precision on z (used to stop the Newton-Raphson method), 
   *   a reasonable choice may be |p2.z - p1.z| / 1000.
   *   Internally, the value can't be higher than |z2 - z1| / 50  nor lower than 1e-15.
   * @param nIterMax upper limit on the number of iteration to be used in the Newton-Raphson method 
   *   a reasonable choice may be 20.
   * @return null or the point such that the tangent line to great-circle arc (if it exists) 
   *   has a slope equals to `+-1`.
   */
  static CooXYZ arcSpecialPointInEqr(final CooXYZ p1, final CooXYZ p2, 
                                            double zEpsMax, final int nIterMax) {
    Vect3D coneCenter =  CooXYZ.crossProd(p1, p2).normalized();
    final double z0 = coneCenter.z();
    final double z1 = p1.z();
    final double z2 = p2.z();
    final boolean north_point = z0 < 0.0;
    // Compute constants
    final double cte = (north_point ? -1.5 : 1.5) * PI_OVER_FOUR;
    // Remark: r = 1 - 2 sin^2(pi/2 / 2) = 1 - 2 * (sqrt(2)/2)^2 = 0
    final double w0 = 1.0 - pow2(z0);
    // Test if we start the method or not.
    final double d1 = fEqr(z1, z0, w0, cte, 0.0);
    final double d2 = fEqr(z2, z0, w0, cte, 0.0);
    if (haveSameSign(d1, d2)) {
      return null;
    }
    // Newton-Raphson method
    zEpsMax = Math.min(zEpsMax, 0.2e-1 * abs(z2 - z1));
    zEpsMax = Math.max(zEpsMax, 1.0e-15);
    double z = 0.5 * (z1 + z2); // mean of z1 and z2
    double zEps = 1.0;
    int nIter = 0;
    while (nIter < nIterMax && abs(zEps) > zEpsMax) {
      zEps = fOverDfEqr(z, z0, w0, cte, 0.0);
      z -= zEps;
      nIter += 1;
    }
    // z must be in the [z1, z2] range
    assert((z1 <= z2 && z1 <= z && z <= z2) || (z2 < z1 &&  z2 <= z && z <= z1));
    if (abs(z) < TRANSITION_Z) {
      final Vect3D v = intersectSmallCircle(p1, p2, z);
      return new CooXYZ(v);
    } else {
      return null;
    }
  }
  
  private static double fEqr(double z, double z0, double w0, double cte, double r) {
    final double w = 1.0 - pow2(z); // in equatortial region, -2/3 < z < 2/3
    final double q = z / w;         // so q is always defined
    final double n = r - z * z0;
    return (z0 - q * n) / sqrt(w0 * w - pow2(n)) - cte;
  }

  private static double fOverDfEqr(double z, double z0, double w0, double cte, double r) {
    final double w = 1 - z * z;
    final double q = z / w;
    final double n = r - z * z0;
    final double sqrtD2MinusN2 = sqrt(w0 * w - n * n);
    final double qn = q * n;
    final double dalphadz = (z0 - qn) / sqrtD2MinusN2;
    final double f = dalphadz - cte;
    final double df = ( q * (2 * z0 - 3 * qn) - n * (1 / w + dalphadz * dalphadz))
        / sqrtD2MinusN2;
    return f / df;
  }

  /**
   * Should not be used: just for debug of the small radius approximation. 
   * @param latStart
   * @param coneCenterLonModHalfPi
   * @param coneCenterLat
   * @param cosOfConeCenterLat
   * @param squareOfSinOfHalfRadius
   * @param positiveSlope
   * @param latEpsMax
   * @param nIterMax
   * @param c
   * @return
   * @throws IOException
   */
  public static double newtonSolveEquatorialZone(final double latStart,
      final double coneCenterLat, final double cosOfConeCenterLat, final double squareOfSinOfHalfRadius,
      final boolean positiveSlope, final double latEpsMax, final int nIterMax, final AngularDistanceComputer c) {
    final double cte = 1.5 * (positiveSlope ? PI_OVER_FOUR : -PI_OVER_FOUR); // 1.5 = 1 / TRANSITION_Z
    double prevLat = 4, lat = latStart;
    int nIter = 0;
    for (int i = 0; i < nIterMax && Math.abs(lat - prevLat) > latEpsMax; i++, nIter++) {
      final double deltaLat = lat - coneCenterLat;
      final double sinDeltaLat = c.sin(deltaLat);
      final double sinHalfDeltaLat = c.sin(0.5 * deltaLat);
      final double cosLat = cos(lat);
      final double tanLat = tan(lat);
      final double n = squareOfSinOfHalfRadius - sinHalfDeltaLat * sinHalfDeltaLat;
      final double d = cosOfConeCenterLat * cosLat;
      final double oneOverd2 = 1 / sqrt(n * (d - n));
      final double oneOverCosLat = 1 / cosLat;
      final double n2 = tanLat * n - 0.5 * sinDeltaLat;
      final double np2 = (n * oneOverCosLat - 0.5 * cosOfConeCenterLat) * oneOverCosLat;
      final double dp2 = 0.5 * (sinDeltaLat * (n - 0.5 * d) - tanLat * d * n) * oneOverd2;
      final double dalpha = n2 * oneOverd2;
      final double d2alpha = np2 * oneOverd2 - dp2 * n2 * oneOverd2 * oneOverd2;
      final double f =  dalpha - cte * cosLat; // * oneOverCosLat
      final double df = tanLat * dalpha + d2alpha; // * oneOverCosLat
      prevLat = lat;
      lat -= f / df;
    }
    //    LOG.debug("Number of iterations: {}", nIter);
    return lat;

  }

  
  
  /**
   * Returns the coordinate z (=sin(lat)) of the point on the cone such that d(DeltaX)/dY = +-1,
   * i.e. the point on the cone having the equation Y = (+-)X + b with the larger (smallest)
   * possible y-intercept (b).
   * The function uses the Newton-Raphson method to solve (DeltaX)/dY -+ 1 = 0.
   * @param zStart the starting point of the Newton-Raphson method. A good choice is
   *        3/2 * sin(lat - 0.5 * rRad) for positive slopes, and
   *        3/2 * sin(lat + 0.5 * rRad) for negative slopes
   * @param coneCenterLonModHalfPi the longitude of the center of the cone, modulo pi/2, in radians
   * @param sinOfConeCenterLat the sine of the latitude of the center of the cone (= z_0)
   * @param coneRadiusRad the raidus of the cone (in radians)
   * @param positiveSlope {@code true} if we look for lines of equations Y = X + b 
   *                      (i.e. south part of the cone),
   *                      {@code false} if we look for lines of equations Y = -X + b
   *                      (e.e. north part of the cone)
   * @param zEpsMax the  
   * @param nIterMax the maximum number of iteration in the Newton-Raphson algorithm
   * 
   * east + positive slope => NE
   * east + negative slope => SE
   * !east + positive slope => SW
   * !east + negative slope => NW
   * 
   * @return
   */
  public static double newtonSolveNorthPolarCapZone(double zStart, boolean eastValue,
      double coneCenterLonModHalfPi,
      double sinOfConeCenterLat, double twoSineOfHalfConeRadius, boolean positiveSlope,
      double zEpsMax, int nIterMax) { // yEpsMax = 10e-8 => angle of ~= 2 mas
    final double cte = (positiveSlope ? 0.5 : -0.5) * PI_OVER_FOUR;
    final double w0 = 1 - sinOfConeCenterLat * sinOfConeCenterLat;
    final double r = (1 - 0.5 * twoSineOfHalfConeRadius * twoSineOfHalfConeRadius);
    final double z0 = sinOfConeCenterLat;
    final double direction = eastValue ? 1 : -1;
    double zEps = 1, z = zStart;
    zEpsMax = Math.max(1e-15, Math.min(zEpsMax, (zStart - sinOfConeCenterLat) * 1e-4));
    for (int i = 0; i < nIterMax && Math.abs(zEps) > zEpsMax; i++) {
      zEps = fOverDfNpc(z, coneCenterLonModHalfPi, z0, w0, cte, direction, r);
      z -= zEps;
    }
    return z;
  }

  
  private static CooXYZ arcSpecialPointInPc(CooXYZ p1, CooXYZ p2, double zEpsMax, int nIterMax) {
    // Ensure p1.lon() < p2.lon()
    if (p1.lon() > p2.lon()) {
      CooXYZ tmp = p1;
      p1 = p2;
      p2 = tmp;
    }
    // Check if great-circle arc overlap several base cells
    assert(p1.lon() % HALF_PI >= 0.0);
    assert(p2.lon() % HALF_PI >= 0.0);
    final int lon1DivHalfPi = (int) (p1.lon() / HALF_PI);
    final int lon2DivHalfPi = (int) (p2.lon() / HALF_PI);
    assert(lon1DivHalfPi < 4);
    assert(lon2DivHalfPi < 4);
    assert(lon1DivHalfPi <= lon2DivHalfPi);
    if (lon1DivHalfPi != lon2DivHalfPi) {
      CooXYZ resZ1 = null;
      CooXYZ resZ2 = null;
      // Info to understand the code
      // - the normal to the plane ra = q * pi/2, with q = 0 or 2 is n = (0, 1, 0)
      // - the normal to the plane ra = q * pi/2, with q = 1 or 3 is n = (1, 0, 0)
      if ((p2.lon() - p1.lon()) > Math.PI) {
        // Cross lon = 0
        Vect3D p2p1N = crossProd(p2, p1).normalized();
        if (p2.lon() % HALF_PI > 0) { // First quarter, [p2.lon, ((p2.lon % PI/2) + 1) * PI/2]
          assert(lon2DivHalfPi > 0);
          int n2Y = lon2DivHalfPi & 1;
          final Vect3D n2 = new Vect3D(n2Y ^ 1, n2Y, 0.0);
          final CooXYZ intersect2 = new CooXYZ(intersectPointPc(p1, p2, p2p1N, n2));
          /*if (intersect2.lon() <= 0.0) {
            assert lon2DivHalfPi == 3;
          }*/
          assert p2.lon() < intersect2.lon() || lon2DivHalfPi == 3: p2.lon() + " " + intersect2.lon();
          resZ2 = arcSpecialPointInPcSameQuarter(p2, intersect2, zEpsMax, nIterMax);
        }
        // Second quarter [(p1.lon % PI/2) * PI/2, p1.lon]
        assert(lon1DivHalfPi < 3);
        int n1X = lon1DivHalfPi & 1;
        final Vect3D n1 = new Vect3D(n1X , n1X ^ 1, 0.0);
        CooXYZ intersect1 = new CooXYZ(intersectPointPc(p1, p2, p2p1N, n1));
        assert intersect1.lon() < p1.lon() : "p1.lon: " + p1.lon() + " vs " + intersect1.lon() + " :int.lon and p2.lon: " + p2.lon();
        resZ1 = arcSpecialPointInPcSameQuarter(intersect1, p1, zEpsMax, nIterMax);
      } else {
        Vect3D p1p2N = crossProd(p1, p2).normalized();
        if (p1.lon() % HALF_PI > 0) { // First quarter, [p1.lon, ((p1.lon % PI/2) + 1) * PI/2]
          assert(lon1DivHalfPi < 3);
          int n1Y = lon1DivHalfPi & 1;
          final Vect3D n1 = new Vect3D(n1Y ^ 1, n1Y, 0.0);
          final CooXYZ intersect1 = new CooXYZ(intersectPointPc(p1, p2, p1p2N, n1));
          assert(p1.lon() < intersect1.lon());
          resZ1 = arcSpecialPointInPcSameQuarter(p1, intersect1, zEpsMax, nIterMax);
        }
        // Last quarter [(p2.lon % PI/2) * PI/2, p2.lon]
        assert(lon2DivHalfPi > 0);
        int n2X = lon2DivHalfPi & 1;
        final Vect3D n2 = new Vect3D(n2X , n2X ^ 1, 0.0);
        final CooXYZ intersect2 = new CooXYZ(intersectPointPc(p1, p2, p1p2N, n2));
        assert(intersect2.lon() < p2.lon());
        resZ2 = arcSpecialPointInPcSameQuarter(intersect2,  p2, zEpsMax, nIterMax);
      }
      if (resZ1 != null) {
        return resZ1;
      } else {
        return resZ2;
      }
    } else { // Same quarter
      return arcSpecialPointInPcSameQuarter(p1, p2, zEpsMax, nIterMax);
    }
  }

// Here we assume that the great-circle arc is in a same quarter
// (i.e. (p1.lon % pi/2) == (p2.lon % pi/2) (except if one of the two point is on a border n * pi/2)  
private static CooXYZ arcSpecialPointInPcSameQuarter(CooXYZ p1, CooXYZ p2, double zEpsMax, int nIterMax) {
  //Ensure p1.lon() < p2.lon()
  if (p1.lon() > p2.lon()) {
    CooXYZ tmp = p1;
    p1 = p2;
    p2 = tmp;
  }
  assert p1.lon() <= p2.lon() : p1.lon() + " " + p2.lon();
  final int p1DivHalfPi = (int) (p1.lon() / HALF_PI);
  double p2ModHalfPi = p2.lon() % HALF_PI;
  if (p2ModHalfPi == 0) {
    p2ModHalfPi = HALF_PI;
  }
  final CooXYZ v1 = new CooXYZ(p1.lon() % HALF_PI, p1.lat());
  final CooXYZ v2 = new CooXYZ(p2ModHalfPi, p2.lat());
  Vect3D coneCenter = crossProd(v1, v2).normalized();
  if (coneCenter.toLon() > PI) {
    coneCenter = coneCenter.opposite();
  }
  
  double coneCenterLon = coneCenter.toLon();
  //debug_assert!(0 <= cone_center_lon && cone_center_lon <= );
  double z0 = coneCenter.z();
  double z1 = v1.z();
  double z2 = v2.z();
  boolean northValue = z0 < 0.0;
  // ( here we could have but do not use the fact that we previously ensure that p1.lon() < p2.lon() )
  final boolean eastValue = ((v1.lat() > v2.lat()) ^ (v1.lon() > v2.lon())) ^ !northValue;
  // Deal with NPC / SPC
  double z = 0.5 * (z1 + z2); // (0.1 * z1 + 0.9 * z2); //(z1 + z2).half(); // mean of z1 and z2
  final boolean spc =  z < 0.0;   // south polar cap
  if (spc) {
    z = -z;
    z1 = -z1;
    z2 = -z2;
    z0 = -z0;
    northValue = !northValue;
  }
  // Compute constants
  //  - remark: r = 1 - 2 sin^2(pi/2 / 2) = 1 - 2 * (sqrt(2)/2)^2 = 0
  final double cte = 0.5 * (northValue ? -PI_OVER_FOUR : PI_OVER_FOUR);
  double w0 = 1 - pow2(z0);
  final double direction = eastValue ? 1.0 : -1.0;
  // Test if we start the method or not
  final double d1 = fNpc(z1, coneCenterLon, z0, w0, cte, direction, 0.0);
  final double d2 = fNpc(z2, coneCenterLon, z0, w0, cte, direction, 0.0);
  if (haveSameSign(d1, d2)) {
    return null;
  }
  // Choose an initial value
  double dz = fOverDfNpc(z, coneCenterLon, z0, w0, cte, direction, 0.0);
  z -= dz;
  if (!((z1 < z && z < z2) || (z2 < z && z < z1))) {
    z = z2 - fOverDfNpc(z2, coneCenterLon, z0, w0, cte, direction, 0.0);
    if (!((z1 < z && z < z2) || (z2 < z && z < z1))) {
      z = z1 - fOverDfNpc(z1, coneCenterLon, z0, w0, cte, direction, 0.0);
    }
  }
  // Newton-Raphson method
  zEpsMax = Math.max(Math.min(zEpsMax, 0.2e-1 * abs(z2 - z1)), 1.0e-15);
  int nIter = 0;
  double zEps = 1.0;
  while (nIter < nIterMax && abs(zEps) > zEpsMax) {
    zEps = fOverDfNpc(z, coneCenterLon, z0, w0, cte, direction, 0.0);
    z -= zEps;
    nIter += 1;
  }
  // Return result if seems correct
  if (/*Double.isFinite(z)*/isFinite(z) && z > TRANSITION_Z && ((z1 < z && z < z2) || (z2 < z && z < z1))) {
    if (spc) {
      final Vect3D v = intersectSmallCircle(p1, p2, -z);
      return new CooXYZ(v);
    } else {
      final Vect3D v = intersectSmallCircle(p1, p2, z);
      return new CooXYZ(v);
    }
  } else {
    return null;
  }
}

/// Returns the intersection point between the given arc (of given normal vector)
/// and the plane of given normal vector
/// WARNING: only valid in polar caps since we use 'z' to determine if we have to take the 
/// result of (p1 x p2) x n or its complements (here we have the guarantee that 
/// sign(p1.z) = sign(p2.z) must be = sign(res.z)
  private static Vect3D intersectPointPc(final CooXYZ p1, final CooXYZ p2,
      final Vect3D p1Xp2, final Vect3D n) {
  assert(abs(p1.z()) >= TRANSITION_Z  && abs(p2.z()) >= TRANSITION_Z);
  assert(p1.z() > 0.0 == p2.z() > 0.0);
  final Vect3D intersect = Vect3D.crossProd(p1Xp2, n).normalized();
  if (!haveSameSign(intersect.z(), p1.z())) {
    return intersect.opposite();
  } else {
    return intersect;
  }
}
  
  
  private static double fNpc(double z, double coneCenterLonModHalfPi, 
      double z0, double w0, double cte, double direction, double r) {
    final double w = 1 - z ;
    final double w2 = 1 - pow2(z);
    final double q = z / w2;
    final double n = r - z * z0;
    final double d2 = w0 * w2;
    final double sqrtD2MinusN2 = Math.sqrt(d2 - n * n);
    final double qn = q * n;
    final double arccos = Math.acos(n / Math.sqrt(d2));
    final double dalphadz = (z0 - qn) / sqrtD2MinusN2;
    return direction * w * dalphadz
        - 0.5 * (direction * arccos + coneCenterLonModHalfPi - PI_OVER_FOUR)
        + cte;
  }

  private static double fOverDfNpc(double z, double coneCenterLonModHalfPi, 
      double z0, double w0, double cte, double direction, double r) {
    final double w = 1 - z ;
    final double w2 = 1 - pow2(z);
    final double q = z / w2;
    final double n = r - z * z0;
    final double d2 = w0 * w2;
    final double sqrtD2MinusN2 = Math.sqrt(d2 - n * n);
    final double qn = q * n;
    final double arccos = Math.acos(n / Math.sqrt(d2));
    final double dalphadz = (z0 - qn) / sqrtD2MinusN2;
    final double f = direction * w * dalphadz
        - 0.5 * (direction * arccos + coneCenterLonModHalfPi - PI_OVER_FOUR)
        + cte;
    final double df = -1.5 * direction * dalphadz + (direction * w / sqrtD2MinusN2) *
        ( q * (2 * z0 - 3 * qn) - n * (1 / w2 + dalphadz * dalphadz));
    return f  /df;
  }

 

  /**
   * Returns the coordinate z (=sin(lat)) of the point on the cone such that d(DeltaX)/dY = +-1,
   * i.e. the point on the cone having the equation Y = (+-)X + b with the larger (smallest)
   * possible y-intercept (b).
   * The function uses the Newton-Raphson method to solve (DeltaX)/dY -+ 1 = 0.
   * @param zStart the starting point of the Newton-Raphson method. A good choice is
   *        3/2 * sin(lat - 0.5 * rRad) for positive slopes, and
   *        3/2 * sin(lat + 0.5 * rRad) for negative slopes
   * @param coneCenterLonModHalfPi the longitude of the center of the cone, modulo pi/2, in radians
   * @param sinOfConeCenterLat the sine of the latitude of the center of the cone (= z_0)
   * @param coneRadiusRad the raidus of the cone (in radians)
   * @param positiveSlope {@code true} if we look for lines of equations Y = X + b 
   *                      (i.e. south part of the cone),
   *                      {@code false} if we look for lines of equations Y = -X + b
   *                      (e.e. north part of the cone)
   * @param zEpsMax the  
   * @param nIterMax the maximum number of iteration in the Newton-Raphson algorithm
   * 
   * east + positive slope => NE
   * east + negative slope => SE
   * !east + positive slope => SW
   * !east + negative slope => NW
   * 
   * @return
   * @throws IOException 
   */
  public static double newtonSolveNorthPolarCapZone(final double latStart, final boolean eastValue,
      final double coneCenterLonModHalfPi, final double coneCenterLat,
      final double cosOfConeCenterLat, final double squareOfSinOfHalfRadius,
      final boolean positiveSlope,
      final double latEpsMax, final int nIterMax, final AngularDistanceComputer c) { //throws IOException {

    /*LOG.debug("START FILE " + String.format(Locale.US, "%.6f", Math.toDegrees(coneCenterLat)));
    writeDAlphaFileNorthPolarCapZon(eastValue,
        coneCenterLonModHalfPi, coneCenterLat,
        cosOfConeCenterLat, squareOfSinOfHalfRadius,
        positiveSlope,
        1000, c, String.format(Locale.US, "%.6f", Math.toDegrees(coneCenterLat)));
     */

    final double cte = positiveSlope ? HALF_PI : 0; // : 0;
    final double direction = eastValue ? 1 : -1;
    // System.err.println("positive slope: " + positiveSlope + "; east direction: " + eastValue);

    double prevLat = 4, lat = latStart;

    int nIter = 0;
    for (int i = 0; i < nIterMax && Math.abs(lat - prevLat) > latEpsMax; i++, nIter++) {
      final double deltaLat = lat - coneCenterLat;
      final double sinDeltaLat = c.sin(deltaLat);
      final double sinHalfDeltaLat = c.sin(0.5 * deltaLat);
      final double cosLat = cos(lat);
      final double tanLat = tan(lat);
      final double n = squareOfSinOfHalfRadius - sinHalfDeltaLat * sinHalfDeltaLat;
      final double d = cosOfConeCenterLat * cosLat;

      final double oneOverd2 = 1 / sqrt(n * (d - n));
      final double oneOverCosLat = 1 / cosLat;
      final double oneOverCosLatminTanLat =  oneOverCosLat - tanLat;
      final double n2 = tanLat * n - 0.5 * sinDeltaLat;
      // final double onePlusTanLatSquare = (1 + tanLat * tanLat);
      final double np2 = (n * oneOverCosLat - 0.5 * cosOfConeCenterLat) * oneOverCosLat;
      final double dp2 = 0.5 * (sinDeltaLat * (n - 0.5 * d) - tanLat * d * n) * oneOverd2;
      final double alpha = coneCenterLonModHalfPi + direction * 2 * c.asin(sqrt(n / d));
      final double dalpha = direction * n2 * oneOverd2;

      final double d2alpha = direction * (np2 * oneOverd2 - dp2 * n2 * oneOverd2 * oneOverd2);

      final double f =  2 * oneOverCosLatminTanLat * dalpha - alpha + cte;
      final double df = 2 * oneOverCosLatminTanLat * (d2alpha - dalpha * oneOverCosLat) - dalpha;

      prevLat = lat;
      lat -= f / df;

      // final double np = -0.5 * sinDeltaLat;
      //final double dp = -d * tanLat; // = -cosOfConeCenterLat * sinLat

      // double check = ((np * d - dp * n) / (d * d)) / (sqrt(n / d) * sqrt(1 - n / d));
      // assert dalpha == check : dalpha + " != " + check; // OK!
      // double check2 = 0.5 * (np * d + n * dp - 2 * n * np) / d2;
      // assert check2 == dp2 : check2 + " != " + dp2;
      // double check3 = tanLat * np + n * oneOverCosLat * oneOverCosLat - 0.5 * cos(deltaLat);
      // assert check3 == np2 : check3 + " != " + np2;

      /*LOG.debug("deltaLat: {}\n"
          + "sinDeltaLat: {}\n"
          + "sinHalfDeltaLat: {}\n"
          + "cosLat: {}\n"
          + "tanLat: {}\n"
          + "n: {}\n"
          + "d: {}\n"
          + "oneOverd2: {}\n"
          // + "d2: {}\n"
          + "oneOverCosLatminTanLat: {}\n"
          + "n2: {}\n"
          //+ "onePlusTanLatSquare: {}\n"(
          + "np2: {}\n"
          + "dp2: {}\n"
          + "alpha: {}\n"
          + "dalpha: {}\n"
          + "d2alpha: {}\n"
          + "f: {}\n"
          + "df: {}\n"
          + "lat: {}\n"
          + "newLat: {}\n",
          deltaLat, 
          sinDeltaLat,
          sinHalfDeltaLat,
          cosLat, 
          tanLat,
          n, 
          d,
          oneOverd2,
          // d2,
          oneOverCosLatminTanLat,
          n2,
          //onePlusTanLatSquare,
          np2,
          dp2,
          alpha,
          dalpha,
          d2alpha,
          f,
          df,
          prevLat,
          lat
          );*/
      // 1 x cos
      // 1 x tan
      // 1 x sin
      // 2 x sqrt
      // 2 x fast sin + 1 x fast cos + 1 x fast arcsin
    }
    //    LOG.debug("Number of iterations (exact): {}", nIter);
    return lat;
  }


  
/// Returns the intersection of the given great-circle arc (defined by the smallest distance 
/// between the two given points) and the small circle of given z (equation $`z=cte`$).
/// Let's use the following notations:
/// - Coordinates of $`\vec{a}\times\vec{b} = (x_0, y_0, z_0)`$
/// - Coordinates of the points we are looking for $`\vec{i} = (x, y, z=cte)`$ 
/// We look for `x` and `y` solving 
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     x^2 + y^2 + z^2 & = & 1 \\
///     xx_0 + yy_0 + zz_0 & = & 0 \\
///     (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 & = & 2 \mathrm(unused)
///   \end{array}
/// \right.
/// ```
/// It leads to
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     y & = & - \left(x\frac{x_0}{y_0} - \frac{zz_0}{y_0}\right) \\
///     0 & = & x^2(1+\frac{x_0^2}{y_0^2}) + x\frac{2x_0zz_0}{y_0^2} + (\frac{zz_0}{y_0})^2 - (x_0^2 + y_0^2)
///   \end{array}
/// \right.
/// ```
/// We solve the quadratic equation
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     ax^2 + bx + c & = & 0 \\
///     \Delta & = & b^2 - 4ac \\
///     x & = & \frac{-b\pm \sqrt{\Delta}}{2a} 
///   \end{array}
/// \right.
/// ```
/// If $`y_0 = 0`$, we directly derive the result:
/// ```math
/// \left\{
///   \begin{array}{rcl}
///     x & = & -\frac{zz_0}{x_0} \\
///     y & = & \pm\sqrt{1 - x^2 - z^2} = \pm\sqrt{1 - z^2(1 + \frac{z_0^2}{x_0^2})}
///   \end{array}
/// \right.
/// ```
/// In both cases, two solutions are available.
/// We select the pair $`(x, y)`$ such that both scalar products with the two great-circle arc 
/// points are higher than the two points scalar product.
public static Vect3D intersectSmallCircle(CooXYZ p1, CooXYZ p2, double z) {
  assert(-1.0 <= z && z <= 1.0);
  if  ((p1.z() < z && z < p2.z()) || (p2.z() < z && z < p1.z())) {
    final double p1p2 = p1.scalarProd(p2);
    final Vect3D p1Xp2 = crossProd(p1, p2).normalized();
    final double x0 = p1Xp2.x();
    final double y0 = p1Xp2.y();
    final double z0 = p1Xp2.z();
    if (abs(y0) <= 1e-14) {
      final double  x = -(z * z0) / x0;
      final double  y1 = sqrt(1 - (pow2(x) + pow2(z)));
      final double  y2 = -y1;
      if (p1.x() * x + p1.y() * y1 + p1.z() * z >= p1p2
        && p2.x() * x + p2.y() * y1 + p2.z() * z >= p1p2) {
        return new Vect3D(x, y1, z);
      } else if (p1.x() * x + p1.y() * y2 + p1.z() * z >= p1p2
             && p2.x() * x + p2.y() * y2 + p2.z() * z >= p1p2) {
        return new Vect3D(x, y2, z);
      } else {
        assert(false);
        return null;
      }
    } else {
      final double  x0_y0 = x0 / y0;
      final double  zz0_y0 = z * z0 / y0;
      final double  a = 1.0 + pow2(x0_y0);
      final double  b = 2.0 * x0_y0 * zz0_y0;
      final double  c = pow2(zz0_y0) + pow2(z) - 1;
      final double  sqrt_delta = sqrt(pow2(b) - 4.0 * a * c);
      final double  x1 = (-b + sqrt_delta) / twice(a);
      final double  y1 = -x1 * x0_y0 - zz0_y0;
      final double  x2 = (-b - sqrt_delta) / twice(a);
      final double  y2 = -x2 * x0_y0 - zz0_y0;
      if  (p1.x() * x1 + p1.y() * y1 + p1.z() * z >= p1p2
        && p2.x() * x1 + p2.y() * y1 + p2.z() * z >= p1p2) {
        return new Vect3D(x1, y1, z);
      } else if (p1.x() * x2 + p1.y() * y2 + p1.z() * z >= p1p2
              && p2.x() * x2 + p2.y() * y2 + p2.z() * z >= p1p2) {
        return new Vect3D(x2, y2, z);
      } else {
        assert(false);
        return null;
      }
    }
  } else {
    return null;
  }
}

  /// Returns the intersection of the given great-circle arc (defined by the smallest distance 
  /// between the two given points) and the small circle of equation $`z=2/3`$.
  /// (Internally, we simply call [intersect_small_circle](fn.intersect_small_circle.html) with
  /// z = 2/3).
  private static Vect3D intersectWithTransitionLatNpc(final CooXYZ p1, final CooXYZ p2) {
    return intersectSmallCircle(p1, p2, TRANSITION_Z);
  }
  
  /// Returns the intersection of the given great-circle arc (defined by the smallest distance 
  /// between the two given points) and the small circle of equation $`z=-2/3`$.
  /// (Internally, we simply call [intersect_small_circle](fn.intersect_small_circle.html) with
  /// z = -2/3).
  private static Vect3D intersectWithTransitionLatSpc(CooXYZ p1, CooXYZ p2) {
    return intersectSmallCircle(p1, p2, -TRANSITION_Z);
  }

 
  private static boolean haveSameSign(double d1, double d2) {
    // Try extracting sign bits?
    return d1 == 0.0 || d2 == 0.0 || ((d1 > 0.0) == (d2 > 0.0));
  }
  
  private static final double twice(double x) {
    return 2 * x;
  }

  private static final double pow2(double x) {
    return x * x;
  }

  // COPIED To STAY COMPATIBLE WITH JAVA 6
  public static boolean isFinite(double d) {
    return Math.abs(d) <= Double.MAX_VALUE;
  }





  /*
  private static void writeDAlphaFileNorthPolarCapZon(final boolean eastValue,
      final double coneCenterLonModHalfPi, final double coneCenterLat,
      final double cosOfConeCenterLat, final double squareOfSinOfHalfRadius,
      final boolean positiveSlope,
      final int nSteps, final AngularDistanceComputer c, String prefix) throws IOException {
    final List<String> lines = new ArrayList<String>();
    lines.add("lat,dDeltaLon,ddDeltaLon");

    final double coneRadiusRad =  2 * Math.asin(sqrt(squareOfSinOfHalfRadius));

    final double latMin = coneCenterLat - coneRadiusRad;
    final double latMax = coneCenterLat + coneRadiusRad;
    final double step = (latMax - latMin) / (nSteps - 1);

    final double cte = positiveSlope ? HALF_PI : 0;
    final double direction = eastValue ? 1 : -1;

    for (int i = 0; i < nSteps; i++) {
      final double lat = latMin + step * i;

      final double deltaLat = lat - coneCenterLat;
      final double sinDeltaLat = c.sin(deltaLat);
      final double sinHalfDeltaLat = c.sin(0.5 * deltaLat);
      final double cosLat = cos(lat);
      final double tanLat = tan(lat);
      final double n = squareOfSinOfHalfRadius - sinHalfDeltaLat * sinHalfDeltaLat;
      final double d = cosOfConeCenterLat * cosLat;

      final double np = -0.5 * sinDeltaLat;
      final double dp = -d * tanLat; // = -cosOfConeCenterLat * sinLat

      final double d2 = sqrt(n * (d - n));
      final double oneOverCosLat = 1 / cosLat;
      final double oneOverCosLatminTanLat = oneOverCosLat - tanLat;
      final double n2 = tanLat * n - 0.5 * sinDeltaLat;
      // final double onePlusTanLatSquare = (1 + tanLat * tanLat);
      final double np2 = ((n / cosLat) - (0.5 * coneCenterLat)) / cosLat;
      final double dp2 = 0.5 * (np * d + n * dp - 2 * n * np) / d2; //0.5 * (sinDeltaLat * (n - 0.5 * d) - tanLat * d * n) / d2;
      final double alpha = coneCenterLonModHalfPi + direction * 2 * c.asin(sqrt(n / d));
      final double dalpha = direction * n2 / d2;

      //double check = ((np * d - dp * n) / (d * d)) / (sqrt(n / d) * sqrt(1 - n / d));
      //assert dalpha == check : dalpha + " != " + check; // OK!

      final double d2alpha = direction * ((np2 / d2) - ((dp2 * n2) / (d2 * d2)));

      final double f =  2 * oneOverCosLatminTanLat * dalpha - alpha + cte;
      final double df = 2 * oneOverCosLatminTanLat * (d2alpha - dalpha / cosLat) - dalpha;

      lines.add(lat + "," + f + "," + df);

    }
    Files.write(new File(prefix + ".deltaLon." + eastValue + "." + positiveSlope + ".csv").toPath(), lines, StandardOpenOption.CREATE);
  }

  /* TODO: A METTRE DANS UNE CLASSE DE TEST */
  /*private static void writeDAlphaFile(double sinOfConeCenterLat, double coneRadiusRad, int nSteps)
      throws IOException {
    final List<String> lines = new ArrayList<String>();
    lines.add("lat,dDeltaLon,ddDeltaLon");

    final double b = 1.5; // 3/2
    final double b2 = 2.25; // (3/2)^2
    final double a = 4 / Math.PI;

    final double centerLat = Math.asin(sinOfConeCenterLat);
    final double yMin = b * Math.sin(centerLat - coneRadiusRad);
    final double yMax = b * Math.sin(centerLat + coneRadiusRad);
    final double step = (yMax - yMin) / (nSteps - 1);

    final double z0 = sinOfConeCenterLat;

    final double sqrtW0 = Math.sqrt(1 - z0 * z0);

    final double u0 = (1 - 0.5 * coneRadiusRad * coneRadiusRad) / sqrtW0;
    final double v0 = Math.abs(z0) / sqrtW0;


    for (int i = 0; i < nSteps; i++) {
      final double y = yMin + step * i;

      final double sqrt1mW =  Math.abs(y) / b;
      final double w = 1 - sqrt1mW * sqrt1mW;
      final double p = v0 - u0 * sqrt1mW;
      final double m = u0 - v0 * sqrt1mW;
      final double o = Math.sqrt(w - m * m);
      final double f = Math.copySign(p / (w * o), y) - a;
      final double df = (2 * sqrt1mW * p) / (w * w * o) - u0 / (w * o)
          + (p / (w * o * o * o)) * (sqrt1mW - v0 * m);

      final double lat = Math.asin(y / 1.5);
      lines.add(lat + "," + f + "," + df);

    }
    Files.write(new File("dDeltaLon.text.csv").toPath(), lines, StandardOpenOption.CREATE);
  }*/
  /* 
   private static void writeDAlphaFile(double sinOfConeCenterLat, double coneRadiusRad, 
       boolean positiveSlope, int nSteps) throws IOException {
     final List<String> lines = new ArrayList<String>();
     lines.add("lat,dDeltaLon,ddDeltaLon");

     final double a = 1.5 / (positiveSlope ? FOUR_OVER_PI : - FOUR_OVER_PI);
     final double w0 = 1 - sinOfConeCenterLat * sinOfConeCenterLat;
     final double twoSineOfHalfConeRadius = 2 * Math.sin(0.5 * coneRadiusRad);
     final double r = (1 - 0.5 * twoSineOfHalfConeRadius * twoSineOfHalfConeRadius);
     final double z0 = sinOfConeCenterLat;

     final double centerLat = Math.asin(sinOfConeCenterLat);
     final double zMin = Math.sin(centerLat - coneRadiusRad);
     final double zMax = Math.sin(centerLat + coneRadiusRad);
     final double step = (zMax - zMin) / (nSteps - 1);

     for (int i = 0; i < nSteps; i++) {
       final double z = zMin + step * i;

       final double oneMinusZ2 = 1 -z * z;
       final double q = z / oneMinusZ2;
       final double n = r - z * z0;
       final double d2 = w0 * oneMinusZ2;
       final double f = (z0 - q * n) / Math.sqrt(d2 - n * n) - a;
       final double df = 0;

       final double lat = Math.asin(z);
       lines.add(lat + "," + f + "," + df);

     }
     Files.write(new File("dDeltaLon.text.v2" + positiveSlope + ".csv").toPath(), lines, StandardOpenOption.CREATE);
   }
   */   
  /**
   * 
   * @param args
   * @throws IOException 
   */
  /*public static void main(final String[] args) throws IOException {
    final double lon = Math.toRadians(32); // 098.46467 +40.18486
    final double lat = Math.toRadians(40.18486);

    final double rRad = Math.toRadians(52.49 / 60);

    System.out.println("latMin: " + Math.toDegrees(lat - rRad));
    System.out.println("latMax: " + Math.toDegrees(lat + rRad));


    final double yStartMax = 1.5 * Math.sin(lat + rRad * 0.9);
    final double yStartMin = 1.5 * Math.sin(lat - rRad * 0.9);
    System.out.println("latStart: " + Math.toDegrees(Math.asin(yStartMax / 1.5)));

    long l1 = System.nanoTime();
    final double yMax = newtonSolveEquatorialZone(yStartMax, Math.sin(lat), rRad, true, rRad / 1000, 100);
    final double yMin = newtonSolveEquatorialZone(yStartMin, Math.sin(lat), rRad, false, rRad / 1000, 100);
    long l2 = System.nanoTime();
    System.out.println("latMaxDeg: " + Math.toDegrees(Math.asin(yMax / 1.5))
      + "; latMinDeg: " + Math.toDegrees(Math.asin(yMin / 1.5))
      + "; in " + (l2 - l1) * 1e-6 + " ms" );

    // writeDAlphaFile(Math.sin(lat), rRad, 1000);

    writeDAlphaFile(Math.sin(lat), rRad, true, 1000);
    writeDAlphaFile(Math.sin(lat), rRad, false, 1000);


    final double twoSineOfHalfConeRadius = 2 * Math.sin(0.5 * rRad);
    final double zStartMin = Math.sin(lat - rRad * 0.9);
    final double zStartMax = Math.sin(lat + rRad * 0.9);
    l1 = System.nanoTime();
    final double zMax = newtonSolveEquatorialZone(zStartMax, Math.sin(lat), twoSineOfHalfConeRadius, false, rRad / 1000, 15);
    System.out.println("---");
    final double zMin = newtonSolveEquatorialZone(zStartMin, Math.sin(lat), twoSineOfHalfConeRadius, true, rRad / 1000, 15);
    l2 = System.nanoTime();
    System.out.println("latMaxDeg: " + Math.toDegrees(Math.asin(zMax))
      + "; latMinDeg: " + Math.toDegrees(Math.asin(zMin))
      + "; in " + (l2 - l1) * 1e-6 + " ms" );


  }*/


}
