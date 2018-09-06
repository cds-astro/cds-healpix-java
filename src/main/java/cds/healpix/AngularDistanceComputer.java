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

import cds.healpix.common.math.Math;

import static java.lang.Math.PI;
import static java.lang.Math.sqrt;

import cds.healpix.common.math.TaylorSeries;

abstract class AngularDistanceComputer {

  /**
   * Possible use of a Taylor serie.
   * @param smallAngRad angle not larger than the typical distance this object will have to compute.
   * @return the cosine of the given small angle.
   */
  protected abstract double cos(double smallAngRad);
  /**
   * Possible use of a Taylor serie.
   * @param smallAngRad angle not larger than the typical distance this object will have to compute.
   * @return the sine of the given small angle
   */
  protected abstract double sin(double smallAngRad);
  /**
   * Possible use of a Taylor serie.
   * @param smallValue value of the sine of an angle not larger than the typical distance
   *        this object will have to compute.
   * @return the arcsine, in radians, of the given small value 
   */
  protected abstract double asin(double smallValue);
  
  
  public final double haversineDistInRad(final double deltaLonRad, final double deltaLatRad,
      final double cosLatA, final double cosLatB) {
    return 2 * asin(sqrt(
        squareOfsinOfhalfDistInRad2(deltaLonRad, deltaLatRad, cosLatA, cosLatB)));
  }
  
  public final double squareOfsinOfhalfDistInRad2(final double deltaLonRad, final double deltaLatRad,
      final double cosLatA, final double cosLatB) {
    final double d = sin(0.5 * deltaLatRad);
    final double a = sin(0.5 * deltaLonRad);
    return d * d + a * a * cosLatA * cosLatB;
  }

  /**
   * Return a value between 0 and PI, or NaN if deltaLatRad is larger than the radius of the cone.
   * We recall that any comparison (<, >, <=, >=, ==) other than (!=)  with NaN returns false 
   * @param squareOfSinOfHalfRadius
   * @param cosOfConeCentreLat
   * @param deltaLatRad
   * @return
   */
  public final double coneDeltaLon(final double squareOfSinOfHalfRadius,
      final double cosOfConeCentreLat, final double deltaLatRad, final double sinLat) {
    final double d = sqrt(1 - sinLat * sinLat) * cosOfConeCentreLat;
    double n = sin(0.5 * deltaLatRad);
    n = squareOfSinOfHalfRadius - n * n;
    if (n < 0) {
      return Double.NaN;
    } else if (n >= d) {
      return PI;
    } else {
      return 2 * asin(sqrt(n / d)); // return NaN if n or d is NaN
    }
  }

  /**
   * The implementation of the returned object depends on the typical angle it will have to
   * deal with. 
   * @param angleScaleInRad approximative size of Delta lon, Delta lat and angular distance
   * @return an object computing angular distance
   */
  public static final AngularDistanceComputer getComputer(final double angleScaleInRad) {
    assert angleScaleInRad >= 0;
    // 4.8 10^-06 rad = 1 arcsec
    // 4.8 10^-09 rad = 1 mas 
    // 4.8 10^-12 rad = 1 micro arcsec 
    // 4.8 10^-15 rad = 1 nano  arcsec
    if (angleScaleInRad <= 1e-5) { // ~ 2 arcsec
      return new AngularDistanceComputer() {
        @Override
        protected double cos(double x) {
          // assert !Double.isFinite(x) || x < 1e-4 : x;
          return x < 1e-4 ? TaylorSeries.cosO2(x) : Math.cos(x);
        }
        @Override
        protected double sin(double x) {
          // assert !Double.isFinite(x) || x < 1e-4 : x;
          return x < 1e-4 ? TaylorSeries.sinO3(x) : Math.sin(x);
        }
        @Override
        protected double asin(double x) {
          // assert !Double.isFinite(x) || x < 1e-4 : x;
          return x < 1e-4 ? TaylorSeries.asinO3(x) : Math.asin(x);
        }
      };
    } else if (angleScaleInRad <= 1e-3) { // ~ 3.4 arcmin 
      return new AngularDistanceComputer() {
        @Override
        protected double cos(double x) {
          // assert !Double.isFinite(x) || x < 1e-2 : x;
          return x < 1e-2 ?TaylorSeries.cosO4(x) : Math.cos(x);
        }
        @Override
        protected double sin(double x) {
          // assert !Double.isFinite(x) || x < 1e-2 : x;
          return x < 1e-2 ? TaylorSeries.sinO5(x) : Math.sin(x);
        }
        @Override
        protected double asin(double x) {
          // assert !Double.isFinite(x) || x < 1e-2 : x;
          return x < 1e-2 ? TaylorSeries.asinO5(x) : Math.asin(x);
        }
      };
    } else if (angleScaleInRad <= 1e-2) { // ~ 34 arcmin 
      return new AngularDistanceComputer() {
        @Override
        protected double cos(double x) {
          // assert !Double.isFinite(x) || x < 0.5e-1: x;
          return x < 0.5e-1 ? TaylorSeries.cosO6(x) : Math.cos(x);
        }
        @Override
        protected double sin(double x) {
          // assert !Double.isFinite(x) || x < 0.5e-1 : x;
          return x < 0.5e-1 ? TaylorSeries.sinO7(x) : Math.sin(x);
        }
        @Override
        protected double asin(double x) {
          // assert !Double.isFinite(x) || x < 0.5e-1 : x;
          return x < 0.5e-1 ? TaylorSeries.asinO7(x) : Math.asin(x);
        }
      };
    } else {
      return new AngularDistanceComputer() {
        @Override
        protected double cos(double x) {
          return Math.cos(x);
        }
        @Override
        protected double sin(double x) {
          return Math.sin(x);
        }
        @Override
        protected double asin(double x) {
          return Math.asin(x);
        }
      };
    }
  }
  
}
