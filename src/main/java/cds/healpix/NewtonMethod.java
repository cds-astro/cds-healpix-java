// Copyright 2017-2018 - Université de Strasbourg/CNRS
// The CDS HEALPix library is developped by the Centre de Données
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

import static cds.healpix.Healpix.TRANSITION_Z;
import static cds.healpix.common.math.Math.FOUR_OVER_PI;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.sqrt;
import static cds.healpix.common.math.Math.tan;

import java.io.File;
import java.io.IOException;
// import java.nio.file.Files;
// import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

// import org.slf4j.Logger;
// import org.slf4j.LoggerFactory;

final class NewtonMethod {
  
//  private static final Logger LOG = LoggerFactory.getLogger("fr.unistra.cds.healpix.cone");

  private static final double SQRT_2_OVER_3 = 0.81649658092772603272;
  
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
      final double o = Math.sqrt(w - m * m);
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
    final double a = 1.5 * (positiveSlope ? PI_OVER_FOUR : -PI_OVER_FOUR); // 1.5 = 1 / TRANSITION_Z
    final double w0 = 1 - sinOfConeCenterLat * sinOfConeCenterLat;
    final double r = (1 - 0.5 * twoSineOfHalfConeRadius * twoSineOfHalfConeRadius);
    final double z0 = sinOfConeCenterLat;
    
    double prevZ = -1, z = zStart;
    int nIter = 0;
    for (int i = 0; i < nIterMax && Math.abs(z - prevZ) > zEpsMax; i++, nIter++) {
      final double oneMinusZ2 = 1 - z * z;
      final double q = z / oneMinusZ2;
      final double n = r - z * z0;
      final double d2 = w0 * oneMinusZ2;
      final double sqrtD2MinusN2 = Math.sqrt(d2 - n * n);
      final double qn = q * n;
      final double dalphadz = (z0 - qn) / sqrtD2MinusN2;
      final double f = dalphadz - a;
      final double df = ( q * (2 * z0 - 3 * qn) - n * (1 / oneMinusZ2 + dalphadz * dalphadz))
          / sqrtD2MinusN2;
      prevZ = z;
      z -= f / df;
    }
//    LOG.debug("Number of iterations: {}", nIter);
    return z;
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
    
    double prevZ = -1, z = zStart;

    zEpsMax = Math.max(1e-15, Math.min(zEpsMax, (zStart - sinOfConeCenterLat) * 1e-4));

    int nIter = 0;
    for (int i = 0; i < nIterMax && Math.abs(z - prevZ) > zEpsMax; i++, nIter++) {
      final double oneMinusZ = 1 - z ;
      final double oneMinusZ2 = 1 - z * z;
      final double q = z / oneMinusZ2;
      final double n = r - z * z0;
      final double d2 = w0 * oneMinusZ2;
      final double sqrtD2MinusN2 = Math.sqrt(d2 - n * n);
      final double qn = q * n;
      final double arccos = Math.acos(n / Math.sqrt(d2));
      final double dalphadz = (z0 - qn) / sqrtD2MinusN2;
      final double f = direction * oneMinusZ * dalphadz
          - 0.5 * (direction * arccos + coneCenterLonModHalfPi - PI_OVER_FOUR)
          + cte;
      final double df = -1.5 * direction * dalphadz + (direction * oneMinusZ / sqrtD2MinusN2) *
          ( q * (2 * z0 - 3 * qn) - n * (1 / oneMinusZ2 + dalphadz * dalphadz));
      prevZ = z;
/* LOG.debug("z: {}; oneMinusZ: {}, oneMinusZ2: {}; q: {}; n: {}; d2: {}; sqrtD2MinusN2: {}; qn: {}; arccos: {}; dalphadz: {}, f: {}; df: {}; zEpsMax: {}",
    z, oneMinusZ, oneMinusZ2, q, n, d2, sqrtD2MinusN2, qn, arccos, dalphadz, f, df, zEpsMax); */
      z -= f / df;
    }
//    LOG.debug("Number of iterations: {}", nIter);
    return z;
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
