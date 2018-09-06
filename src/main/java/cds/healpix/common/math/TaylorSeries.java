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
 * Utility class continaing the Taylor expansion of a few trigonometric functions.
 *
 * We recall that:
 *   arcos(x) = pi/2 - arcsin(x)
 * 
 * 
 * @author F.-X. Pineau
 *
 */
public final class TaylorSeries {

  /** Constant value equal to 1/2. */
  public static final double HALF = 0.5d;

  /** Constant value equal to 1/3. */
  public static final double ONE_THIRD = 0.33333333333333333333d; // 1 / 3d;

  /** Constant value equal to 1/6 = 1/6. */
  public static final double ONE_SIXTH = 0.16666666666666666666d; // 1 / 6d;

  /** Constant value equal to 1/24 = 1/4!. */
  public static final double ONE_OVER_24 = 0.04166666666666666666; // 1 / 24d;

  /** Constant value equal to 1/120 = 1/5!. */
  public static final double ONE_OVER_120 = 0.00833333333333333333d; //1 / 120d;

  /** Constant value equal to 1/720 = 1/6!. */
  public static final double ONE_OVER_720 = 0.00138888888888888888d; // 1 / 720d;

  /** Constant value equal to 1/5040 = 1/7!. */
  public static final double ONE_OVER_5040 = 0.00019841269841269841d; // 1 / 5040d;

  /** Constant value equal to 2/15. */
  public static final double TWO_OVER_15 = 0.13333333333333333333d;// 2 / 15d;

  /** Constant value equal to 3/40. */
  public static final double THREE_OVER_40 = 0.075d; // 3 / 40d;

  /** Constant value equal to 3/40. */
  public static final double FIVE_OVER_112 = 0.04464285714285714285d; // 5 / 112d;

  /** Constant value equal to 17/215. */
  public static final double Q_17_OVER_315 = 0.05396825396825396825d; // 17 / 315d;

  /**
   * Returns the Taylor serie of the trigonometric cosine around a = 0 at the 0 order (constant).
   * That is: {@code 1}.
   * The error on this approximation is no more than {@code |x|^2 / 2},  i.e. ~5e-17 if x = 1e-8.
   * We recall that the erro due to numercial approximation is 1e-16.
   * @param a angle, in radians
   * @return 1, i.e. the Taylor serie of the trigonometric cosine around a = 0 at the 0 order.
   */
  public static final double cosO0(final double a) {
    // assert a <= 1e-8;
    return 1;
  }

  /**
   * Returns the Taylor serie of the trigonometric cosine around a = 0 at the order of the 2nd degree.
   * That is: {@code 1 - x^2/2}.
   * The error on this approximation is no more than {@code |x|^4 / 4!},  i.e. ~4e-18 if x = 1e-4.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric sine around a = 0 at the order of the 2nd degree.
   */
  public static final double cosO2(final double a) {
    // assert a <= 1e-4;
    return 1 - HALF * a * a; // + O(a^4) (max 4e-18 if a <= 1e-4)
  }

  /**
   * Returns the Taylor serie of the trigonometric cosine around a = 0 at the order of the 4th degree.
   * That is: {@code 1 - x^2/2 + x^4/4!}.
   * The error on this approximation is no more than {@code |x|^6 / 6!},  i.e. ~1.3e-21 if x = 1e-3.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric sine around a = 0 at the order of the 4th degree.
   */
  public static final double cosO4(final double a) {
    // assert a <= 1e-3;
    final double a2 = a * a;
    return 1 - HALF * a2 + ONE_OVER_24 * a2 * a2; // + O(a^6) (max 1.3e-21 if a <= 1e-3)
  }

  /**
   * Returns the Taylor serie of the trigonometric cosine around a = 0 at the order of the 6th degree.
   * That is: {@code 1 - x^2/2 + x^4/4! - x^6/6!}.
   * The error on this approximation is no more than {@code |x|^8 / 8!},  i.e. ~1.7e-20 if x = 1e-2.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric sine around a = 0 at the order of the 6th degree.
   */
  public static final double cosO6(final double a) {
    // assert a <= 1e-2;
    final double a2 = a * a;
    final double a4 = a2 * a2;
    return 1 - HALF * a2 + ONE_OVER_24 * a4 - ONE_OVER_720 * a2  * a4; // + O(a^8) (max 1.7e-20 if a <= 1e-2)
  }



  /**
   * Returns the Taylor serie of the trigonometric sine around a = 0 at the order of the 1st degree.
   * That is: {@code x}.
   * The error on this approximation is no more than {@code |x|^3 / 3!},  i.e. ~1.6e-19 if x = 1e-6.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric sine around a = 0 at the order of the 1st degree.
   */
  public static final double sinO1(final double a) {
    // assert a <= 1e-6;
    return a;
  }

  /**
   * Returns the Taylor serie of the trigonometric sine around a = 0 at the order of the 3rd degree.
   * That is: {@code x - x^3/3!}.
   * The error on this approximation is no more than {@code |x|^5 / 120}, i.e. ~8.3e-18 if x = 1e-3.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric sine around a = 0 at the order the 3rd degree.
   */
  public static final double sinO3(final double a) {
    // assert a <= 1e-4;
    // For a better precision, we do not use x * (1 + x * x / 6), idem for all other Taylor series
    return a - ONE_SIXTH * a * a * a; // + O(a^5) (max 1e-25 if a <= 1e-5)
  }

  /**
   * Returns the Taylor serie of the trigonometric sine around a = 0 at the order of the 5th degree.
   * That is: {@code x - x^3/3! + x^5/5!}.
   * The error on this approximation is no more than {@code |x|^7 / 7!}, i.e. ~2e-18 if x = 1e-2.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric sine around a = 0
   * at the order of the 5th degree.
   */
  public static final double sinO5(final double a) {
    // assert a <= 1e-3;
    final double a2 = a * a;
    final double a3 = a * a2;
    return a - ONE_SIXTH * a3 + ONE_OVER_120 * a2 * a3; // + O(a^7) (max 2e-24 if a <= 1e-3)
  }

  /**
   * Returns the Taylor serie of the trigonometric sine around a = 0 at the order of the 7th degree.
   * That is: {@code x - x^3/3! + x^5/5! - x^7/7!}.
   * The error on this approximation is no more than {@code |x|^9 / 9!}, i.e. ~2.8e-15 if x = 1e-1.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric sine around a = 0
   * at the order of the 7th degree.
   */
  public static final double sinO7(final double a) {
    // assert a <= 1e-2;
    final double a2 = a * a;
    final double a3 = a * a2;
    final double a5 = a2 * a3;
    return a - ONE_SIXTH * a3 + ONE_OVER_120 * a5 - ONE_OVER_5040 * a2 * a5; // + O(a^9) (max 3e-24 if a <= 1e-2)
  }

  /**
   * Returns the Taylor serie of the trigonometric tan around a = 0 at the order of the 1st degree.
   * That is: {@code x}.
   * The error on this approximation is no more than {@code |x|^3 / 3}, i.e. ~4e-19 if x = 1e-6.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric tan around a = 0 at the order of the 1st degree.
   */
  public static final double tanO1(final double a) {
    // assert a <= 1e-6;
    return a;
  }

  /**
   * Returns the Taylor serie of the trigonometric tan around a = 0 at the order of the 3rd degree.
   * That is: {@code x + x^3/3}.
   * The error on this approximation is no more than {@code 2*|x|^5 / 15}, i.e. ~1.3e-21 if x = 1e-4.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric tan around a = 0 at the order of the 3rd degree.
   */
  public static final double tanO3(final double a) {
    // assert a <= 1e-4;
    return a + ONE_THIRD * a * a * a;
  } 

  /**
   * Returns the Taylor serie of the trigonometric tan around a = 0 at the order of the 5th degree.
   * That is: {@code x + x^3/3 + 2x^5/15}.
   * The error on this approximation is no more than {@code 17|x|^7 / 315}, i.e. ~5.4e-23 if x = 1e-3.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric tan around a = 0 at the order of the 5th degree.
   */
  public static final double tanO5(final double a) {
    // assert a <= 1e-3;
    final double a2 = a * a;
    final double a3 = a * a2;
    return a + ONE_THIRD * a3 + TWO_OVER_15 * a2 * a3;
  } 

  /**
   * Returns the Taylor serie of the trigonometric tan around a = 0 at the order of the 7th degree.
   * That is: {@code x + x^3/3 + 2x^5/15 + 17x^7/315}.
   * The error on this approximation is no more than {@code 62|x|^9 / 2835}, i.e. ~2.2e-20 if x = 1e-2.
   * @param a angle, in radians
   * @return the Taylor serie of the trigonometric tan around a = 0 at the order of the 7th degree.
   */
  public static final double tanO7(final double a) {
    // assert a <= 1e-2;
    final double a2 = a * a;
    final double a3 = a * a2;
    final double a5 = a2 * a3;
    return a + ONE_THIRD * a3 + TWO_OVER_15 * a5 + Q_17_OVER_315 * a2 * a5;
  } 

  public static final double asinO1(final double x) {
    return x;
  }

  public static final double asinO3(final double x) {
    final double x2 = x * x;
    return x + ONE_SIXTH * x * x2;
  }

  public static final double asinO5(final double x) {
    final double x2 = x * x;
    final double x3 = x * x2;
    return x + ONE_SIXTH * x3 + THREE_OVER_40 * x2 * x3;
  }

  public static final double asinO7(final double x) {
    final double x2 = x * x;
    final double x3 = x * x2;
    final double x5 = x2 * x3;
    return x + ONE_SIXTH * x3 + THREE_OVER_40 * x5 + FIVE_OVER_112 * x2 * x5;
  }

  // atan

}
