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

package cds.healpix.common.math;

/**
 * Utility class to easily switch between java.lang.Math and another Math
 * librarie like FastMath, and to define new functions like sinc
 * (sinus cardinal = sin(x) / x).
 * 
 * Replace java..lang.Math by org..apache.commons.math3.util.FastMath
 * (using e.g. sed or replaceAll) to use FastMath instead of Math
 * (and vice-versa).
 *
 * @author F.-X. Pineau
 *
 */
public final class Math {

  public static final double EPSILON = 1e-15;
  
  // Constants in 1 / n;
  public static final double ONE_THIRD = 1d / 3d; // 1 / 3
  public static final double ONE_FOURTH = 0.25d;  // 1 / 4
  
  // Constants including PI.
  public static final double PI = java.lang.Math.PI;
  // - in PI / n
  public static final double HALF_PI = 0.5 * PI;             // PI / 2
  public static final double PI_OVER_FOUR = ONE_FOURTH * PI; // PI / 4
  // - in n * PI
  public static final double TWO_PI = 2 * PI;  // 2 PI
  // - in n / PI
  public static final double ONE_OVER_PI = 1 / PI;      // 1 / PI
  public static final double ONE_OVER_2PI = 0.5 * ONE_OVER_PI; 
  public static final double TWO_OVER_PI = 2 * ONE_OVER_PI; // 2 / PI
  public static final double FOUR_OVER_PI = 4 * ONE_OVER_PI;

  public static final double LOG2 = log(2);
  
  public static final double SQRT2 = sqrt(2);
  public static final double SQRT2_INV = 1 / SQRT2;
  public static final double SQRT3 = sqrt(3);
  public static final double SQRT3_INV = 1 / SQRT3;
  public static final double SQRT6 = sqrt(6);


  public static boolean isFinite(final double v){
    return Math.abs(v) <= Double.MAX_VALUE;
  }

  public static double toRadians(final double a) {
      return java.lang.Math.toRadians(a);
  }

  public static double toDegrees(final double a) {
      return java.lang.Math.toDegrees(a);
  }

  public static double abs(final double a) {
      return java.lang.Math.abs(a);
  }
  
  public static double signum(final double a) {
      return java.lang.Math.signum(a);
  }

  public static int round(final float a) {
      return java.lang.Math.round(a);
  }
  
  public static long round(final double a) {
      return java.lang.Math.round(a);
  }
  
  public static double pow(final double a, double b) {
      return java.lang.Math.pow(a, b);
  }

  public static double cos(final double a) {
      return java.lang.Math.cos(a);
  }

  public static double sin(final double a) {
      return java.lang.Math.sin(a);
  }

  /*public static double oneMinusSin(double a) {
    a *= 0.5;
    
  }*/
  
  /**
   * To be used when sin(a) is near from 1 and a high precision is needed on 1 - sin(a).
   * Uses formula:
   *   2 * cos(x / 2 + pi/4)
   * Coming from:
   *      1 - sin(x)
   *   =  sin(pi / 2) - sin(x)
   *   =  2 sin(pi / 4 - x / 2) cos(pi / 4 + x / 2)
   *   = -2 sin(x / 2 - pi / 4) cos(x / 2 + pi / 4)
   *   = -2 sin(x / 2 + pi / 4 - pi / 2) cos(x / 2 + pi / 4)
   *   =  2 cos(x / 2 + pi / 4) cos(x / 2 + pi / 4)
   * @param a angle in radians
   * @return 1 - sin(a)
   */
  public static double oneMinusSin(double a) {
    a *= 0.5;
    a += PI_OVER_FOUR;
    a = cos(a);
    a *= a;
    a *= 2;
    return a; // 2 * cos^2(a/2 + pi/4)
    // Should I write in this more readable form?
    // a = cos(0.5 * a - PI_OVER_FOUR);
    // return 2 * a * a;
  }
  
  /**
   * To be used when sin(a) is near from 1 and a high precision is needed on sqrt(1 - sin(a)).
   * {@link #oneMinusSin}.
   * @param a angle in radians HAVING a POSITIVE COSINE.
   * @return sqrt(1 - sin(a))
   */
  public static double sqrtOfOneMinusSin(double a) {
    a *= 0.5;
    a += PI_OVER_FOUR;
    a = cos(a);
    return SQRT2 * (a > 0 ? a : -a);
  }
  
  /**
   * Same as {@link #sqrtOfOneMinusSin} but limited to angle having a positive cosine
   * (in practice cos(a/2 + PI/2) must be positive).
   * @param a angle in radians HAVING a POSITIVE COSINE, or more generally, a in [-3pi/2, pi/2].
   * @return sqrt(1 - sin(a))
   */
  public static double sqrtOfOneMinusSinPC(double a) {
    assert -3 * HALF_PI <= a && a <= HALF_PI;
    a *= 0.5;
    a += PI_OVER_FOUR;
    a = cos(a);
    a *= SQRT2;
    return a;
  }
  
  public static double tan(final double a) {
      return java.lang.Math.tan(a);
  }

  /**
   * Returns the cardinal sine function of the given angle, i.e. sin(x) / x.
   * Precision of 10e-16 on small angles.
   * @param a angle in radians
   * @return the cardinal sine function of the given angle, i.e. sin(x) / x.
   */
  public static double sinc(double a) {
      if (abs(a) > 1.e-4) {
          return sin(a) / a;
      } else {
          // If a is mall, use Taylor expension of asin(a) / a
          // a = 1e-4 => a^4 = 1.e-16
          a *= a;
          return 1 - a * (1 - a / 20.0) / 6.0;
      }
  }

  /**
   * Same as {@link #sinc(double)} but assuming the argument is positive.
   * Precision of 10e-16 on small angles.
   * @param a angle we are looking for the sine, must be &gt; 0.
   * @return the cardinal sine function of the given angle,
   * assuming it is positive.
   */
  public static double sincP(double a) {
      assert a >= 0;
      if (a > 1.e-4) {
          return sin(a) / a;
      } else {
          // If a is mall, use Taylor expension of asin(a) / a
          // a = 1e-4 => a^4 = 1.e-16
          a *= a;
          return 1 - a * (1 - a / 20.0) / 6.0;
      }
  }

  public static double acos(final double a) {
      return java.lang.Math.acos(a);
  }

  public static double asin(final double a) {
      return java.lang.Math.asin(a);
  }

  /**
   * Returns the inverse of the cardinal sine function of the given argument,
   * i.e. the inverse of sin(x) / x.
   * @param a argument
   * @return the inverse of the cardinal sine function of the given aegument.
   */
  public static double asinc(double a) {
      if (abs(a) > 1.e-4) {
          return asin(a) / a;
      } else {
          // If a is mall, use Taylor expension of asin(a) / a
          // a = 1e-4 => a^4 = 1.e-16
          a *= a;
          return 1 + a * (1 + a * 9.0 / 20.0) / 6.0; 
      }
  }

  /**
   * Same as {@link #asinc(double)} but assuming the argument is positive.
   * @param a argument
   * @return the inverse of the cardinal sine function of the given aegument,
   * assuming the argument is positive.
   */
  public static double asincP(double a) {
      assert a >= 0;
      if (a > 1.e-4) {
          return asin(a) / a;
      } else {
          // If a is mall, use Taylor expension of asin(a) / a
          // a = 1e-4 => a^4 = 1.e-16
          a *= a;
          return 1 + a * (1 + a * 9.0 / 20.0) / 6.0;
      }
  }

  public static double atan2(final double y, final double x) {
      return java.lang.Math.atan2(y, x);
  }

  public static double tanh(final double x) {
      return java.lang.Math.tan(x);
  }

  /**
   * Returns tanh-1(x), i.e. the inverse function of tanh.
   * @param x argument, in range ]-1, 1[ (NaN returned otherwise)
   * @return the hyperbolic inverse tangent of the given argument
   */
  public static double atanh(final double x) {
      return 0.5 * log((1 + x) / (1 - x));
  }

  public static double log(final double x) {
      return java.lang.Math.log(x);
  }

  public static double sqrt(final double angleRadians) {
      return java.lang.Math.sqrt(angleRadians);
  }

  /**
   * Compute the number of real roots of the quadratic equation:
   * q2*x^2 + q1*x + q0 = 0
   * @param q0 coefficient of x^0
   * @param q1 coefficient of x^1
   * @param q2 coefficient of x^2
   * @param result arrays of size at least 2 storing the real roots (if any)
   * @return the number of real roots
   */
  public static int roots(final double q0, final double q1, final double q2,
          final double[] result) {
      double delta = q1 * q1 - 4 * q2 * q0; // b^2 - 4ac
      int i = 0;
      if (delta == 0) {
          result[i++] = -((0.5 * q1) / q2) ; // -b / 2a
      } else if (delta > 0) {
          delta = sqrt(delta);
          result[i++] = -0.5 * (delta + q1) / q2; // (-b - sqrt(b^2 - 4ac)) / 2a
          result[i++] =  0.5 * (delta - q1) / q2; // (-b + sqrt(b^2 - 4ac)) / 2a
      } // else assert delta < 0 || Double.isNaN(delta)
      return i;
  }

  public static int roots(double q0, double q1, double q2, double q3,
          final double[] result) {
      // Reduce to x^3 + p x^2 + q x + r = 0
      // p      q         r
      q2 /= q3; q1 /= q3; q0 /= q3;
      // That reduces to y^3 + a y + b = 0 with x = y - p/3
      // Compute a = (3q - p^2) / 3
      q3 = (3 * q1 - q2 * q2) / 3; 
      // Compute b = (2p^3 - 9pq + 27r) / 27
      q2 = (2 * q2 * q2 * q2 - 9 * q1 * q2 + 27 * q0) / 27;
      // Compute b^2/4 + a^3/27
      double delta = (0.25 * q2 * q2) + (q3 * q3 * q3) / 27;
      int i = 0;
      if (delta > 0) {
          delta = sqrt(delta);
          q0 = -0.5 * q2; // -b / 2
          // solution y = A + B
          // with A = [-b/2 - sqrt(b^2/4 + a^3/27)]^1/3
          //      B = [-b/2 + sqrt(b^2/4 + a^3/27)]^1/3
          q0 = pow(q0 - delta, ONE_THIRD) + pow(q0 + delta, ONE_THIRD);
          result[i++] = q0 - ONE_THIRD * q1;
      } else if (delta == 0) {
          if (q2 < 0) {
              q3 = sqrt(-q3 / 3);
              result[i++] = 2 * q3;
              result[i++] = -q3;
              result[i++] = -q3;
          } else if (q2 ==0) {
              result[i++] = 0;
              result[i++] = 0;
              result[i++] = 0;
          } else {
              assert q2 > 0 : q2;
              q3 = sqrt(-q3 / 3);
              result[i++] = -2 * q3;
              result[i++] = q3;
              result[i++] = q3;
          }
      } else if (delta < 0) {
          throw new Error("TO BE CONTINUED!!");
      } else {
          assert !isFinite(delta);
      }
      assert 0 <= i && i <= 3;
      return i;
  }
  
  /*public static double normalizeFromMinusPiToPlusPi(double a) {
    a %= TWO_PI;                // to range ]-2*PI, 2*PI[  
    a = (a + TWO_PI) % TWO_PI;  // to range [0, 2*PI[  
    if (a > PI) {               // to range ]-PI, PI] 
       a -= TWO_PI; 
    }
    return a;
    See code in apache FastMath!!
  }*/

  // lat % 180 => ]-180, 180[, if (>90 ou < -90, ra+= 180 (do not chnage the sinus)
  
}

