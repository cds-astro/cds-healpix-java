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
 * Comes from the Apache FastMath class, and is thus under the Apache licence.
 *
 * @author F.-X. Pineau
 *
 */
public final class FastMath {

  /**
   * 0x40000000 - used to split a double into two parts, both with the low order bits cleared.
   * Equivalent to 2^30.
   */
  private static final long HEX_40000000 = 0x40000000L; // 1073741824L

  /** Mask used to clear low order 30 bits */
  private static final long MASK_30BITS = -1L - (HEX_40000000 -1); // 0xFFFFFFFFC0000000L;
  
  /** Constant: {@value}. */
  private static final double F_1_2 = 1d / 2d;

  /** Eighths.
   * This is used by sinQ, because its faster to do a table lookup than
   * a multiply in this time-critical routine
   */
  private static final double EIGHTHS[] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625};

  /** Sine, Cosine, Tangent tables are for 0, 1/8, 2/8, ... 13/8 = PI/2 approx. */
  private static final int SINE_TABLE_LEN = 14;

  /** Sine table (high bits). */
  private static final double SINE_TABLE_A[] =
    {
        +0.0d,
        +0.1246747374534607d,
        +0.24740394949913025d,
        +0.366272509098053d,
        +0.4794255495071411d,
        +0.5850973129272461d,
        +0.6816387176513672d,
        +0.7675435543060303d,
        +0.8414709568023682d,
        +0.902267575263977d,
        +0.9489846229553223d,
        +0.9808930158615112d,
        +0.9974949359893799d,
        +0.9985313415527344d,
    };

  /** Sine table (low bits). */
  private static final double SINE_TABLE_B[] =
    {
        +0.0d,
        -4.068233003401932E-9d,
        +9.755392680573412E-9d,
        +1.9987994582857286E-8d,
        -1.0902938113007961E-8d,
        -3.9986783938944604E-8d,
        +4.23719669792332E-8d,
        -5.207000323380292E-8d,
        +2.800552834259E-8d,
        +1.883511811213715E-8d,
        -3.5997360512765566E-9d,
        +4.116164446561962E-8d,
        +5.0614674548127384E-8d,
        -1.0129027912496858E-9d,
    };

  /** Cosine table (high bits). */
  private static final double COSINE_TABLE_A[] =
    {
        +1.0d,
        +0.9921976327896118d,
        +0.9689123630523682d,
        +0.9305076599121094d,
        +0.8775825500488281d,
        +0.8109631538391113d,
        +0.7316888570785522d,
        +0.6409968137741089d,
        +0.5403022766113281d,
        +0.4311765432357788d,
        +0.3153223395347595d,
        +0.19454771280288696d,
        +0.07073719799518585d,
        -0.05417713522911072d,
    };

  /** Cosine table (low bits). */
  private static final double COSINE_TABLE_B[] =
    {
        +0.0d,
        +3.4439717236742845E-8d,
        +5.865827662008209E-8d,
        -3.7999795083850525E-8d,
        +1.184154459111628E-8d,
        -3.43338934259355E-8d,
        +1.1795268640216787E-8d,
        +4.438921624363781E-8d,
        +2.925681159240093E-8d,
        -2.6437112632041807E-8d,
        +2.2860509143963117E-8d,
        -4.813899778443457E-9d,
        +3.6725170580355583E-9d,
        +2.0217439756338078E-10d,
    };

  /** Tangent table, used by atan() (high bits). */
  private static final double TANGENT_TABLE_A[] =
      {
      +0.0d,
      +0.1256551444530487d,
      +0.25534194707870483d,
      +0.3936265707015991d,
      +0.5463024377822876d,
      +0.7214844226837158d,
      +0.9315965175628662d,
      +1.1974215507507324d,
      +1.5574076175689697d,
      +2.092571258544922d,
      +3.0095696449279785d,
      +5.041914939880371d,
      +14.101419448852539d,
      -18.430862426757812d,
  };

  /** Tangent table, used by atan() (low bits). */
  private static final double TANGENT_TABLE_B[] =
      {
      +0.0d,
      -7.877917738262007E-9d,
      -2.5857668567479893E-8d,
      +5.2240336371356666E-9d,
      +5.206150291559893E-8d,
      +1.8307188599677033E-8d,
      -5.7618793749770706E-8d,
      +7.848361555046424E-8d,
      +1.0708593250394448E-7d,
      +1.7827257129423813E-8d,
      +2.893485277253286E-8d,
      +3.1660099222737955E-7d,
      +4.983191803254889E-7d,
      -3.356118100840571E-7d,
  };
  
  // Generic helper methods
  
  public static final double SAFE_MIN = 0x1.0p-1022;
  
  /**
   * Get the high order bits from the mantissa.
   * Equivalent to adding and subtracting HEX_40000 but also works for very large numbers
   *
   * @param d the value to split
   * @return the high order part of the mantissa
   */
  private static double doubleHighPart(double d) {
      if (d > -SAFE_MIN && d < SAFE_MIN){
          return d; // These are un-normalised - don't try to convert
      }
      long xl = Double.doubleToRawLongBits(d); // can take raw bits because just gonna convert it back
      xl &= MASK_30BITS; // Drop low order bits
      return Double.longBitsToDouble(xl);
  }

  
  
  /**
   *  Computes sin(x) - x, where |x| &lt; 1/16.
   *  Use a Remez polynomial approximation.
   *  @param x a number smaller than 1/16
   *  @return sin(x) - x
   */
  private static double polySine(final double x)
  {
    double x2 = x*x;

    double p = 2.7553817452272217E-6;
    p = p * x2 + -1.9841269659586505E-4;
    p = p * x2 + 0.008333333333329196;
    p = p * x2 + -0.16666666666666666;
    //p *= x2;
    //p *= x;
    p = p * x2 * x;

    return p;
  }

  /**
   *  Computes cos(x) - 1, where |x| &lt; 1/16.
   *  Use a Remez polynomial approximation.
   *  @param x a number smaller than 1/16
   *  @return cos(x) - 1
   */
  private static double polyCosine(double x) {
    double x2 = x*x;

    double p = 2.479773539153719E-5;
    p = p * x2 + -0.0013888888689039883;
    p = p * x2 + 0.041666666666621166;
    p = p * x2 + -0.49999999999999994;
    p *= x2;

    return p;
  }


  /**
   *  Compute sine over the first quadrant (0 &lt; x &lt; pi/2).
   *  Use combination of table lookup and rational polynomial expansion.
   *  @param xa number from which sine is requested
   *  @param xb second param
   *  @return sin(xa + xb)
   */
  public static double sinQ(double xa, double xb) {
    int idx = (int) ((xa * 8.0) + 0.5);
    final double epsilon = xa - EIGHTHS[idx]; //idx*0.125;

    // Table lookups
    final double sintA = SINE_TABLE_A[idx];
    final double sintB = SINE_TABLE_B[idx];
    final double costA = COSINE_TABLE_A[idx];
    final double costB = COSINE_TABLE_B[idx];

    // Polynomial eval of sin(epsilon), cos(epsilon)
    double sinEpsA = epsilon;
    double sinEpsB = polySine(epsilon);
    final double cosEpsA = 1.0;
    final double cosEpsB = polyCosine(epsilon);

    // Split epsilon   xa + xb = x
    final double temp = sinEpsA * HEX_40000000;
    double temp2 = (sinEpsA + temp) - temp;
    sinEpsB +=  sinEpsA - temp2;
    sinEpsA = temp2;

    /* Compute sin(x) by angle addition formula */
    double result;

    double a = 0;
    double b = 0;

    double t = sintA;
    double c = a + t;
    double d = -(c - a - t);
    a = c;
    b += d;

    t = costA * sinEpsA;
    c = a + t;
    d = -(c - a - t);
    a = c;
    b += d;

    b = b + sintA * cosEpsB + costA * sinEpsB;

    b = b + sintB + costB * sinEpsA + sintB * cosEpsB + costB * sinEpsB;

    if (xb != 0.0) {
      t = ((costA + costB) * (cosEpsA + cosEpsB) -
          (sintA + sintB) * (sinEpsA + sinEpsB)) * xb;  // approximate cosine*xb
      c = a + t;
      d = -(c - a - t);
      a = c;
      b += d;
    }

    result = a + b;

    return result;
  }

  /**
   * Compute cosine in the first quadrant by subtracting input from PI/2 and
   * then calling sinQ.  This is more accurate as the input approaches PI/2.
   *  @param xa number from which cosine is requested
   *  @return cos(xa + xb)
   */
  public static double cosQ(double xa) {
    final double pi2a = 1.5707963267948966;
    final double pi2b = 6.123233995736766E-17;

    final double a = pi2a - xa;
    double b = -(a - pi2a + xa);
    b += pi2b;

    return sinQ(a, b);
  }


  /** Compute the arc sine of a number.
   * @param x number on which evaluation is done
   * @return arc sine of x
   */
  public static double asin(double x) {
    if (x != x) {
      return Double.NaN;
    }

    if (x > 1.0 || x < -1.0) {
      return Double.NaN;
    }

    if (x == 1.0) {
      return Math.PI/2.0;
    }

    if (x == -1.0) {
      return -Math.PI/2.0;
    }

    if (x == 0.0) { // Matches +/- 0.0; return correct sign
      return x;
    }

    /* Compute asin(x) = atan(x/sqrt(1-x*x)) */

    /* Split x */
    double temp = x * HEX_40000000;
    final double xa = x + temp - temp;
    final double xb = x - xa;

    /* Square it */
    double ya = xa*xa;
    double yb = xa*xb*2.0 + xb*xb;

    /* Subtract from 1 */
    ya = -ya;
    yb = -yb;

    double za = 1.0 + ya;
    double zb = -(za - 1.0 - ya);

    temp = za + yb;
    zb += -(temp - za - yb);
    za = temp;

    /* Square root */
    double y;
    y = sqrt(za);
    temp = y * HEX_40000000;
    ya = y + temp - temp;
    yb = y - ya;

    /* Extend precision of sqrt */
    yb += (za - ya*ya - 2*ya*yb - yb*yb) / (2.0*y);

    /* Contribution of zb to sqrt */
    double dx = zb / (2.0*y);

    // Compute ratio r = x/y
    double r = x/y;
    temp = r * HEX_40000000;
    double ra = r + temp - temp;
    double rb = r - ra;

    rb += (x - ra*ya - ra*yb - rb*ya - rb*yb) / y;  // Correct for rounding in division
    rb += -x * dx / y / y;  // Add in effect additional bits of sqrt.

    temp = ra + rb;
    rb = -(temp - ra - rb);
    ra = temp;

    return atan(ra, rb, false);
  }

  /**
   * Arctangent function
   *  @param x a number
   *  @return atan(x)
   */
  public static double atan(double x) {
    return atan(x, 0.0, false);
  }

  /** Internal helper function to compute arctangent.
   * @param xa number from which arctangent is requested
   * @param xb extra bits for x (may be 0.0)
   * @param leftPlane if true, result angle must be put in the left half plane
   * @return atan(xa + xb) (or angle shifted by {@code PI} if leftPlane is true)
   */
  private static double atan(double xa, double xb, boolean leftPlane) {
    if (xa == 0.0) { // Matches +/- 0.0; return correct sign
      return leftPlane ? copySign(Math.PI, xa) : xa;
    }

    final boolean negate;
    if (xa < 0) {
      // negative
      xa = -xa;
      xb = -xb;
      negate = true;
    } else {
      negate = false;
    }

    if (xa > 1.633123935319537E16) { // Very large input
      return (negate ^ leftPlane) ? (-Math.PI * F_1_2) : (Math.PI * F_1_2);
    }

    /* Estimate the closest tabulated arctan value, compute eps = xa-tangentTable */
    final int idx;
    if (xa < 1) {
      idx = (int) (((-1.7168146928204136 * xa * xa + 8.0) * xa) + 0.5);
    } else {
      final double oneOverXa = 1 / xa;
      idx = (int) (-((-1.7168146928204136 * oneOverXa * oneOverXa + 8.0) * oneOverXa) + 13.07);
    }

    final double ttA = TANGENT_TABLE_A[idx];
    final double ttB = TANGENT_TABLE_B[idx];

    double epsA = xa - ttA;
    double epsB = -(epsA - xa + ttA);
    epsB += xb - ttB;

    double temp = epsA + epsB;
    epsB = -(temp - epsA - epsB);
    epsA = temp;

    /* Compute eps = eps / (1.0 + xa*tangent) */
    temp = xa * HEX_40000000;
    double ya = xa + temp - temp;
    double yb = xb + xa - ya;
    xa = ya;
    xb += yb;

    //if (idx > 8 || idx == 0)
    if (idx == 0) {
      /* If the slope of the arctan is gentle enough (< 0.45), this approximation will suffice */
      //double denom = 1.0 / (1.0 + xa*tangentTableA[idx] + xb*tangentTableA[idx] + xa*tangentTableB[idx] + xb*tangentTableB[idx]);
      final double denom = 1d / (1d + (xa + xb) * (ttA + ttB));
      //double denom = 1.0 / (1.0 + xa*tangentTableA[idx]);
      ya = epsA * denom;
      yb = epsB * denom;
    } else {
      double temp2 = xa * ttA;
      double za = 1d + temp2;
      double zb = -(za - 1d - temp2);
      temp2 = xb * ttA + xa * ttB;
      temp = za + temp2;
      zb += -(temp - za - temp2);
      za = temp;

      zb += xb * ttB;
      ya = epsA / za;

      temp = ya * HEX_40000000;
      final double yaa = (ya + temp) - temp;
      final double yab = ya - yaa;

      temp = za * HEX_40000000;
      final double zaa = (za + temp) - temp;
      final double zab = za - zaa;

      /* Correct for rounding in division */
      yb = (epsA - yaa * zaa - yaa * zab - yab * zaa - yab * zab) / za;

      yb += -epsA * zb / za / za;
      yb += epsB / za;
    }

    epsA = ya;
    epsB = yb;

    /* Evaluate polynomial */
    final double epsA2 = epsA * epsA;

    /*
  yb = -0.09001346640161823;
  yb = yb * epsA2 + 0.11110718400605211;
  yb = yb * epsA2 + -0.1428571349122913;
  yb = yb * epsA2 + 0.19999999999273194;
  yb = yb * epsA2 + -0.33333333333333093;
  yb = yb * epsA2 * epsA;
     */

    yb = 0.07490822288864472;
    yb = yb * epsA2 - 0.09088450866185192;
    yb = yb * epsA2 + 0.11111095942313305;
    yb = yb * epsA2 - 0.1428571423679182;
    yb = yb * epsA2 + 0.19999999999923582;
    yb = yb * epsA2 - 0.33333333333333287;
    yb = yb * epsA2 * epsA;


    ya = epsA;

    temp = ya + yb;
    yb = -(temp - ya - yb);
    ya = temp;

    /* Add in effect of epsB.   atan'(x) = 1/(1+x^2) */
    yb += epsB / (1d + epsA * epsA);

    final double eighths = EIGHTHS[idx];


    //result = yb + eighths[idx] + ya;
    double za = eighths + ya;
    double zb = -(za - eighths - ya);
    temp = za + yb;
    zb += -(temp - za - yb);
    za = temp;

    double result = za + zb;

    if (leftPlane) {
      // Result is in the left plane
      final double resultb = -(result - za - zb);
      final double pia = 1.5707963267948966 * 2;
      final double pib = 6.123233995736766E-17 * 2;

      za = pia - result;
      zb = -(za - pia + result);
      zb += pib - resultb;

      result = za + zb;
    }


    if (negate ^ leftPlane) {
      result = -result;
    }

    return result;
  }
  /** Compute the arc cosine of a number.
   * @param x number on which evaluation is done
   * @return arc cosine of x
   */
  public static double acos(double x) {
    if (x != x) {
        return Double.NaN;
    }

    if (x > 1.0 || x < -1.0) {
        return Double.NaN;
    }

    if (x == -1.0) {
        return Math.PI;
    }

    if (x == 1.0) {
        return 0.0;
    }

    if (x == 0) {
        return Math.PI/2.0;
    }

    /* Compute acos(x) = atan(sqrt(1-x*x)/x) */

    /* Split x */
    double temp = x * HEX_40000000;
    final double xa = x + temp - temp;
    final double xb = x - xa;

    /* Square it */
    double ya = xa*xa;
    double yb = xa*xb*2.0 + xb*xb;

    /* Subtract from 1 */
    ya = -ya;
    yb = -yb;

    double za = 1.0 + ya;
    double zb = -(za - 1.0 - ya);

    temp = za + yb;
    zb += -(temp - za - yb);
    za = temp;

    /* Square root */
    double y = sqrt(za);
    temp = y * HEX_40000000;
    ya = y + temp - temp;
    yb = y - ya;

    /* Extend precision of sqrt */
    yb += (za - ya*ya - 2*ya*yb - yb*yb) / (2.0*y);

    /* Contribution of zb to sqrt */
    yb += zb / (2.0*y);
    y = ya+yb;
    yb = -(y - ya - yb);

    // Compute ratio r = y/x
    double r = y/x;

    // Did r overflow?
    if (Double.isInfinite(r)) { // x is effectively zero
        return Math.PI/2; // so return the appropriate value
    }

    double ra = doubleHighPart(r);
    double rb = r - ra;

    rb += (y - ra*xa - ra*xb - rb*xa - rb*xb) / x;  // Correct for rounding in division
    rb += yb / x;  // Add in effect additional bits of sqrt.

    temp = ra + rb;
    rb = -(temp - ra - rb);
    ra = temp;

    return atan(ra, rb, x<0);
  }

  
  
  /** Compute the square root of a number.
   * <p><b>Note:</b> this implementation currently delegates to {@link Math#sqrt}
   * @param a number on which evaluation is done
   * @return square root of a
   */
  public static double sqrt(final double a) {
      return Math.sqrt(a);
  }
  /**
   * Returns the first argument with the sign of the second argument.
   * A NaN {@code sign} argument is treated as positive.
   *
   * @param magnitude the value to return
   * @param sign the sign for the returned value
   * @return the magnitude with the same sign as the {@code sign} argument
   */
  public static double copySign(double magnitude, double sign){
      // The highest order bit is going to be zero if the
      // highest order bit of m and s is the same and one otherwise.
      // So (m^s) will be positive if both m and s have the same sign
      // and negative otherwise.
      final long m = Double.doubleToRawLongBits(magnitude); // don't care about NaN
      final long s = Double.doubleToRawLongBits(sign);
      if ((m^s) >= 0) {
          return magnitude;
      }
      return -magnitude; // flip sign
  }

  
}
