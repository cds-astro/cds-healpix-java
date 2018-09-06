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
 * Utility class which takes its name from the book 
 * "Hacker's Delight" by Henry S. Warren. Jr.
 * We also uses the "Bit Twiddling Hacks" page by Sean Eron Anderson from Standford.
 * 
 * @author F.-X Pineau
 *
 */
public final class HackersDelight {

  /** Bit mask to isolate the sign bit of an int or a float. */
  public static final int SIGN_BIT_MASK_I = 0x80000000;
  
  /** Bit mask to isolate the sign bit of a long or double. */
  public static final long SIGN_BIT_MASK_L     = 0x8000000000000000L;
  /** Bit mask to isolate all but the sign bit of a long or double. */
  public static final long BUT_SIGN_BIT_MASK_L = 0x7FFFFFFFFFFFFFFFL;
  
  /** Bit mask to isolate the sign bit of a long or double. */
  public static final long EXPONENT_BITS_MASK     = 0x7FF0000000000000L;
  /** Bit mask to isolate all but the sign bit of a long or double. */
  public static final long BUT_EXPONENT_BIT_MASK = 0x800FFFFFFFFFFFFFL;
  
  /**
   * Determines if an integer is a power of two, exact solution.
   * From the "Bit Twiddling Hacks" web page of Sean Eron Anderson, slightly modified.
   * @param x the argument
   * @return true if the argument is a power of 2 and is different from 0. 
   */
  public static boolean isPowOf2(int x) {
    return (x != 0) && isPowOf2Fast(x);
  }
  
  /**
   * Determines if an integer is a power of two, including 0.
   * From the "Bit Twiddling Hacks" web page of Sean Eron Anderson, slightly modified.
   * @param x the argument
   * @return if the argument is a power of 2 OR if it is 0 (which is not a power of 2).
   * For an exact solution (excluding 0), {@link #isPowOf2(int)}. 
   */
  public static boolean isPowOf2Fast(int x) {
    return (x & --x) == 0;
  }

  public long flp2(long x) {
    x |= x >>> 1;
    x |= x >>> 2;
    x |= x >>> 4;
    x |= x >>> 8;
    x |= x >>> 16;
    x |= x >>> 32;
    return x - (x >>> 1);
  }
  
  public static double halfOf(final double x) {
    if (x == 0) {
      return x;
    }
    final long bits = toBits(x);
    final long signAndFraction = (bits & BUT_EXPONENT_BIT_MASK);
    long exponent = bits & EXPONENT_BITS_MASK;
    exponent >>= 52;
    --exponent;
    exponent <<= 52;
    return fromBits(signAndFraction | exponent);
  }
  
  public static long toBits(final double x) {
    return Double.doubleToRawLongBits(x);
  }
  
  public static double fromBits(final long doubleBits) {
    return Double.longBitsToDouble(doubleBits);
  }
  
  public static int abs(int x) {
    int mask = x >> 31;
    assert (x < 0 && mask == -1) || (x >= 0 && mask == 0);
    return (x ^ mask) - mask;
  }

  public static long abs(long x) {
    long mask = x >> 63;
    assert (x < 0 && mask == -1L) || (x >= 0 && mask == 0L);
    return (x ^ mask) - mask;
  }
  
  public static double abs(double x) {
    return fromBits(toBits(x) & BUT_SIGN_BIT_MASK_L);
  }
  
  
  public static int floorIntP(final double a) {
    return (int) a;
  }
  
  public static int floorIntN(final double a) {
    int d = (int) a;
    return (d == a) ? d : --d;
  }
  
  public static int floorInt(final double a) { // A TESTER!!
    int d = (int) a;
    return (a < 0 && d != a) ? --d : d;
  }
  
  public static long floorLongP(final double a) {
    return (long) a;
  }
  
  public static long floorLongN(final double a) {
    long d = (long) a;
    return (d == a) ? d : --d;
  }
  
  public static long floorLong(final double a) { // A TESTER!!
    long d = (long) a;
    return (a < 0 && d != a) ? --d : d;
  }
  
  /**
   * Doz stands for "difference or zero" (see p. 37 of HD).
   * It returns {@code 0} if {@code x > y}, else returns {@code x - y},
   * so it can be replace by {@code (x - y) > 0 ? x - y : 0 } 
   * @param x an argument
   * @param y another argument
   * @return {@code 0} if {@code x > y}, else returns {@code x - y}.
   */
  public static int doz(int x, int y) {
    x -= y;
    return x & ~(x >> 31);
  }

  /**
   * Doz stands for "difference or zero" (see p. 37 of HD).
   * It returns {@code 0} if {@code x > y}, else returns {@code x - y},
   * so it can be replace by {@code (x - y) > 0 ? x - y : 0 } 
   * @param x an argument
   * @param y another argument
   * @return {@code 0} if {@code x > y}, else returns {@code x - y}.
   */
  public static long doz(long x, long y) {
    x -= y;
    return x & ~(x >> 31);
  }

  /**
   * Returns the smaller of two values.
   * @param x an argument
   * @param y another argument
   * @return the smaller of two values.
   */
  public static int min(int x, int y) {
    return x - doz(x, y);
  }

  /**
   * Returns the smaller of two values.
   * @param x an argument
   * @param y another argument
   * @return the smaller of two values.
   */
  public static long min(long x, long y) {
    return x - doz(x, y);
  }

  /**
   * Returns the greater of two values.
   * @param x an argument
   * @param y another argument
   * @return the greater of two values.
   */
  public static int max(int x, int y) {
    return y + doz(x, y);
  }

  /**
   * Returns the greater of two values.
   * @param x an argument
   * @param y another argument
   * @return the greater of two values.
   */
  public static long max(long x, long y) {
    return y + doz(x, y);
  }

  /**
   * Branch-free version of {@link java.lang.Integer#numberOfLeadingZeros} (also taken from the HD).
   * See "Number of leading zeros, branch-free binary search" in HD.
   * @param x the argument
   * @return the number of leading zeros in the given argument.
   */
  public static int nlz(int x) {
    // Comments are from HD.
    int y = -(x >> 16);      // If left half of nside is 0,
    int m = (y >> 16) & 16;  // set n ) 16. If left half
    int n = 16 - m;          // is nonzero, set n = 0 and
    x >>>= m;                // shift nside right 16.
                             // Now nside is of the form 0000xxxx.
    y = x - 0x100;           // If positions 8-15 are 0,
    m = (y >> 16) & 8;       // add 8 to n and shift nside left 8.
    n += m;
    x <<= m;
    
    y = x - 0x1000;          // If positions 12-15 are 0,
    m = (y >> 16) & 4;       // add 4 to n and shift nside left 4.
    n += m;
    x <<= m;
    
    y = x - 0x4000;          // If positions 14-15 are 0,
    m = (y >> 16) & 2;       // add 2 to n and shift nside left 2.
    n += m;
    x <<= m;
    
    y = x >> 14;            // Set y = 0, 1, 2 or 3.
    m = y & ~(y >> 1);      // Set m = 0, 1, 2 or 2 resp.
    return n - m + 2;
  }
  
  public static void main(String[] args) {
    double d = -0.0;
    System.out.println(1 >> 1);
    
    System.out.println(d);
    System.out.println(halfOf(-0.0));
    System.out.println("----------------");
    System.out.println(halfOf(0.0));
    System.out.println("----------------");
    System.out.println(halfOf(-2.0));
    System.out.println("----------------");
    System.out.println(halfOf(-458.25));
    
    /*System.out.println("-1L: " + Integer.toBinaryString(-1));
    
    System.out.println(doz(1, 2) + " -- " + doz(2, 1));
    System.out.println(doz(15892, 125) + " -- " + doz(125, 15892));

    int x = 15892, y = 125;
    int d = x - y;
    System.out.println("d: " + d + "; bin: " + Integer.toBinaryString(d));
    System.out.println("d >> 31: " + (d >> 31) + "; bin: " + Integer.toBinaryString(d >> 31));
    System.out.println("~(d >> 31) " + (~(d >> 31)) + "; bin: " + Integer.toBinaryString(~(d >> 31)));

    System.out.println("-----------------");

    y = 15892; x = 125;
    d = x - y;
    System.out.println("d: " + d + "; bin: " + Integer.toBinaryString(d));
    System.out.println("d >> 31: " + (d >> 31) + "; bin: " + Integer.toBinaryString(d >> 31));
    System.out.println("~(d >> 31) " + (~(d >> 31)) + "; bin: " + Integer.toBinaryString(~(d >> 31)));*/
  }

}
