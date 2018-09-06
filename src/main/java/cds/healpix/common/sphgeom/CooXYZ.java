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

package cds.healpix.common.sphgeom;

import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.asin;
import static cds.healpix.common.math.Math.atan2;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.sqrt;

import java.util.Random;

/**
 * Defines a coordinates on the unit-sphere, internally both Euclidean and Spherical coordinates
 * are storef.
 * 
 * @author F.-X. Pineau
 *
 */
public class CooXYZ {

  private final double lon;
  private final double lat;

  private final double x;
  private final double y;
  private final double z;

  public CooXYZ(final double lonRad, final double latRad) {
    this.lon = lonRad;
    this.lat = latRad;
    // equa coo to xyz
    final double cosDec = cos(this.lat);
    this.x = cosDec * cos(this.lon);
    this.y = cosDec * sin(this.lon);
    this.z = sin(this.lat);
  }

  public CooXYZ(final double x, final double y, final double z) {
    this.x = x;
    this.y = y;
    this.z = z;
    // xyz to equa coo
    final double tmp = atan2(this.y, this.x);
    this.lon = tmp < 0 ? 2 * PI + tmp : tmp;
    this.lat = atan2(this.z, sqrt(this.x * this.x + this.y * this.y));
  }

  public static CooXYZ toEquaCooXYZ(final CooXYZ pos) {
    return new CooXYZ(pos.lon(), pos.lat());
  } 

  /**
   * Getter
   * @return the longitude, in radians
   */
  public final double lon() { return this.lon; }
  /**
   * Getter
   * @return the latitude, in radians
   */
  public final double lat() { return this.lat; }

  /**
   * Getter
   * @return the x cartesian coordinate
   */
  public final double x() { return this.x; }
  /**
   * Getter
   * @return the y cartesian coordinate
   */
  public final double y() { return this.y; }
  /**
   * Getter
   * @return the z cartesian coordinate
   */
  public final double z() { return this.z; }

  /**
   * Returns the sum of the given vector, normalized so that the resulting vector is on the unit
   * sphere.
   * @param vects vector we are looking for the normalized sum.
   * @return the sum of the given vector, normalized so that the resulting vector is on the unit
   * sphere.
   */
  public static CooXYZ normalizedSum(final CooXYZ... vects) {
    double x = 0, y = 0, z = 0;
    for (final CooXYZ v : vects) {
      x += v.x;
      y += v.y;
      z += v.z;
    }
    final double norm = sqrt(x * x + y * y + z * z);
    return new CooXYZ(x / norm, y / norm, z / norm);
  }
  
  /**
   * Returns the cross-product of the two given vectors.
   * @param v1 first vector
   * @param v2 second vector
   * @return the cross-product of the two given vectors.
   */
  public static Vect3D crossProd(final CooXYZ v1, final CooXYZ v2) {
    return new Vect3D((v1.y * v2.z) - (v1.z * v2.y)
        ,(v1.z * v2.x) - (v1.x * v2.z)
        ,(v1.x * v2.y) - (v1.y * v2.x));
  }
  
  /**
   * Computes the scalar product of this point with given vectors.
   * @param v the second vector used in the scalar product.
   * @return the scalar product of this point with given vectors.
   */
  public double scalarProd(final Vect3D v) {
    return (x * v.x()) + (y * v.y()) + (z * v.z());
  }
  
  /**
   * Computes the scalar product of this point with given vectors.
   * @param v the second vector used in the scalar product.
   * @return the scalar product of this point with given vectors.
   */
  public double scalarProd(final CooXYZ v) {
    return x * v.x + y * v.y + z * v.z;
  }
  
  /**
   * Returns the spherival distance (using the Haversine formula) separating the two given points.
   * @param c1 first point
   * @param c2 second point
   * @return the spherival distance (using the Haversine formula) separating the two given points.
   */
  public static final double havDist(final CooXYZ c1,
      final CooXYZ c2) {
    final double d = sin(0.5 * (c2.lat() - c1.lat()));
    final double a = sin(0.5 * (c2.lon() - c1.lon()));
    final double h3 = d * d + a * a * cos(c2.lat()) * cos(c1.lat());
    return 2 * asin(sqrt(h3));
  }

  /**
   * Retruns the spherical distance separating the two given points.
   * @param c1 first point
   * @param c2 second point
   * @return the spherical distance separating the two given points.
   */
  public static final double spheDist(final CooXYZ c1,
      final CooXYZ c2) {
    return 2 * asin(0.5 * euclDist(c1, c2));
  }
  /**
   * Returns the Euclidean distance separating the two given points.
   * @param c1 first point
   * @param c2 second point
   * @return the Euclidean distance separating the two given points.
   */
  public static final double euclDist(final CooXYZ c1,
      final CooXYZ c2) {
    double d = 0;
    double tmp = c2.x() - c1.x();
    d += tmp * tmp;
    tmp = c2.y() - c1.y();
    d += tmp * tmp;
    tmp = c2.z() - c1.z();
    d += tmp * tmp;
    return sqrt(d);
  }
  @Override
  public String toString() {
    return "lonDeg: " + Math.toDegrees(lon) + "; latDeg: " + Math.toDegrees(lat) + "; " +  x + "," + y + "," + z;
  }

  /**
   * Returns the minimum enclosing cone, i.e. the cone containig the two given points and having
   * the smallest possible radius. In this trivial case, the diameter of the cone is the arc (ab).
   * @param a first point
   * @param b secobd point
   * @return the minimum enclosing cone, i.e. the cone containig the two given points and having
   * the smallest possible radius. In this trivial case, the diameter of the cone is the arc (ab).
   */
  public static Cone mec(final CooXYZ a, final CooXYZ b) {
    final double r = 0.5 * spheDist(a, b);
    return new Cone(arcCenter(a, b), r);
  }

  /**
   * Returns the Minimum Enclosing Cone, i.e. the cone containig the three given points and having
   * the smallest possible radius.
   * @param a first point
   * @param b secobd point
   * @param c third point
   * @return the Minimum Enclosing Cone, i.e. the cone containig the three given points and having
   * the smallest possible radius.
   */
  public static Cone mec(final CooXYZ a,
      final CooXYZ b, final CooXYZ c) {
    final double da = spheDist(b, c);
    final double db = spheDist(a, c);
    final double dc = spheDist(a, b);
    double r = dc;
    CooXYZ centerPos;
    final boolean only2points = (da > db && da > dc)
        ? (r = 0.5 * da) >= spheDist(centerPos = arcCenter(b, c, r), a)
        : (db > dc) ? (r = 0.5 * db) >= spheDist(centerPos = arcCenter(a, c, r), b)
          : (r = 0.5 * dc) >= spheDist(centerPos = arcCenter(a, b, r), c);
     return only2points ? new Cone(centerPos, r)
         : new Cone(circumCenter(a, b, c, r = circumRadiusSphe(da, db, dc)), r);
  }

  /**
   * Returns the Minimum Enclosing Cone, i.e. the cone containig all the given points and having the
   * smallest possible radius.
   * @param p list of the points we look for the minimum enclising cone
   * @return the Minimum Enclosing Cone, i.e. the cone containig all the given points and having the
   * smallest possible radius.
   */
  public static Cone mec(final CooXYZ... p) {
    Cone cone = null;
    switch(p.length) {
    case 0:
      // cone = new Cone(0, 0, 0);
      break;
    case 1:
      cone = new Cone(p[0], 0);
      break;
    case 2:
      double r = 0.5 * spheDist(p[0], p[1]);
      cone = new Cone(arcCenter(p[0], p[1]), r);
      break;
    case 3:
      // a and b will contains the most distant points
      // c will contains the third point
      CooXYZ a = p[0], b = p[1], c = p[2];
      double da = spheDist(b, c);
      double db = spheDist(a, c);
      double dc = spheDist(a, b);
      if (da > db && da > dc) { // da greatest distance a <--> c
        final CooXYZ ctmp = a;
        final double dtmp = da;
        a = c;
        da = dc;
        c = ctmp;
        dc = dtmp;
      } else if (db > dc) { // db greatest distance b <--> c
        final CooXYZ ctmp = b;
        final double dtmp = db;
        b = c;
        db = dc;
        c = ctmp;
        dc = dtmp;
      }
      r = 0.5 * dc;
      CooXYZ centerPos = arcCenter(a, b, r);
      if (spheDist(c, centerPos) <= r) {
        cone = new Cone(centerPos, r);
      } else {
        r = circumRadiusSphe(da, db, dc);
        cone = new Cone(circumCenter(a, b, c, r) , r);
      }
      break;
    default:
      cone = minSphericalCircle(p.clone());
    }
    return cone;
  }

  /**
   * Returns the angular radius (in radians) of the circumcircle of a
   * spherical triangle of given side lengths a, b and c.
   * @param a first size length (in radians)
   * @param b second size length (in radians)
   * @param c third size length (in radians)
   * @return in the angular radius (in radians) of the circumcircle of a
   * spherical triangle of given side lengths a, b and c.
   */
  public static final double circumRadiusSphe(double a, double b, double c) {
    final double sinas2 = sin(0.5 * a);
    final double sinbs2 = sin(0.5 * b);
    final double sincs2 = sin(0.5 * c);
    // Numerator 
    final double n = sinas2 * sinbs2 * sincs2;
    // Denominator
    final double d = (sinas2 + sinbs2 + sincs2) * (sinas2 + sinbs2 - sincs2)
        * (sinas2 - sinbs2 + sincs2) * (sinbs2 + sincs2 - sinas2); 
    return asin(sqrt((4 * n * n) / d));
  }

  /**
   * Returns the angular radius (in radians) of the circumcircle of a
   * spherical triangle of given vertices a, b and c.
   * @param a first vertex
   * @param b second vertex
   * @param c third vertex
   * @return the angular radius (in radians) of the circumcircle of a
   * spherical triangle of given vertices a, b and c.
   */
  public static final double circumRadiusSphe(final CooXYZ a, final CooXYZ b, final CooXYZ c) {
    final double da = spheDist(b, c);
    final double db = spheDist(a, c);
    final double dc = spheDist(a, b);
    return circumRadiusSphe(da, db, dc);
  }

  /**
   * Returns the center on the unit sphere of the circumcircle of a
   * spherical triangle of given vertices a, b and c.
   * @param a first vertex
   * @param b second vertex
   * @param c third vertex
   * @return the center on the unit sphere of the circumcircle of a
   * spherical triangle of given vertices a, b and c.
   */
  public static final CooXYZ circumCenter(final CooXYZ a, final CooXYZ b, final CooXYZ c) {
    return circumCenter(a, b, c, circumRadiusSphe(a, b, c));
  }

  /**
   * Returns the center on the unit sphere of the circumcircle of radius r of
   * a spherical triangle of given vertices a, b and c.
   * @param a first vertex
   * @param b second vertex
   * @param c third vertex
   * @param r spherical radius of the circumcircle
   * @return the center on the unit sphere of the circumcircle of radius r of
   * a spherical triangle of given vertices a, b and c.
   */
  public static final CooXYZ circumCenter(final CooXYZ a, final CooXYZ b, final CooXYZ c, final double r) {
    final double e = 1 - 0.5 * r * r;
    // Simple cramer resolution of AX = E (here A --> X)
    final double d = a.x() * (b.y() * c.z() - c.y() * b.z())
        - b.x() * (a.y() * c.z() - c.y() * a.z())
        + c.x() * (a.y() * b.z() - b.y() * a.z());
    final double x = e * (b.y() * c.z() - c.y() * b.z())
        - e * (a.y() * c.z() - c.y() * a.z())
        + e * (a.y() * b.z() - b.y() * a.z());
    final double y = a.x() * (e * c.z() - e * b.z())
        - b.x() * (e * c.z() - e * a.z())
        + c.x() * (e * b.z() - e * a.z());
    final double z = a.x() * (b.y() * e - c.y() * e)
        - b.x() * (a.y() * e - c.y() * e)
        + c.x() * (a.y() * e - b.y() * e);
    return new CooXYZ(x / d, y / d, z / d);
  }

  /**
   * Returns the center of the arc define by the smallest distance (on the unit sphere) between
   * the two given points.
   * @param a first point
   * @param b second point
   * @return the center of the arc define by the smallest distance (on the unit sphere) between
   * the two given points.
   */
  public static final CooXYZ arcCenter(final CooXYZ a, final CooXYZ b) {
    return arcCenter(a, b, 0.5 * spheDist(a, b));
  }

  /**
   * Faster version of {@link #arcCenter(CooXYZ, CooXYZ)} when we already know the distance between
   * the two given points.
   * @param a first point
   * @param b second point
   * @param r half the distance between a and b
   * @return the center of the arc define by the smallest distance (on the unit sphere, = 2*r)
   * between the two given points.
   */
  public static final CooXYZ arcCenter(final CooXYZ a, final CooXYZ b, final double r) {
    final double e = 1 - 0.5 * r * r;
    final double e3 = 0; // 3rd contraint: cross product = 0
    final double x13 = a.y() * b.z() - a.z() * b.y(), 
        x23 = a.z() * b.x() - a.x() * b.z(),
        x33 = a.x() * b.y() - a.y() * b.x();
    // Simple cramer resolution of AX = E (here A --> X)
    final double d = a.x() * (b.y() * x33 - x23 * b.z())
        - b.x() * (a.y() * x33   - x23   * a.z())
        +   x13 * (a.y() * b.z() - b.y() * a.z());
    final double x = e * (b.y() * x33 - x23 * b.z())
        - e * (a.y() * x33 - x23 * a.z())
        + e3 * (a.y() * b.z() - b.y() * a.z());
    final double y = a.x() * (e * x33 - e3 * b.z())
        - b.x() * (e * x33 - e3 * a.z())
        + x13 * (e * b.z() - e * a.z());
    final double z = a.x() * (b.y() * e3 - x23 * e)
        - b.x() * (a.y() * e3 - x23   * e)
        + x13   * (a.y() * e  - b.y() * e);
    return new CooXYZ(x / d, y / d, z / d);
  }

  private static Cone minSphericalCircle(final CooXYZ[] p) {
    double r = 0.5 * spheDist(p[0], p[1]);
    Cone c = new Cone(arcCenter(p[0], p[1]), r);
    for (int i = 2; i < p.length; i++) {
      if (!c.contains(p[i])) {
        shuffle(p, 0, i); // try without this and compare performances!
        r = 0.5 * spheDist(p[0], p[i]);
        c =  new Cone(arcCenter(p[0], p[i]), r);
        for (int j = 1; j < i; j++) {
          if (!c.contains(p[j])) {
            r = 0.5 * spheDist(p[j], p[i]);
            c =  new Cone(arcCenter(p[j], p[i]), r);
            for (int k = 0; k < j; k++) {
              if (!c.contains(p[k])) {
                c = mec(p[k], p[j], p[i]);
                r = c.radiusRad();
              }
            }
          }
        }
      }
    }
    return c;
  }

  private volatile static Random r;

  private static final void shuffle(final Object[] a, final int from,
      final int to) {
    Random rnd = r;
    if (rnd == null) {
      synchronized(CooXYZ.class) {
        r = rnd = new Random();
      }
    }
    for (int i = (to - from); i > 1; ) {
      swap(a, rnd.nextInt(i + from), --i + from);
    }
  }

  /**
   * Swaps the two specified elements in the specified array.
   * @param a array in which the swap is performed
   * @param i index of the first element to be swapped
   * @param j index of the second element to be swapped
   */
  private static void swap(final Object[] a, final int i, final int j) {
    final Object tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
  }

}
