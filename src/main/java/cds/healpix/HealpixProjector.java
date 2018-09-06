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


import static cds.healpix.Projection.X_INDEX;
import static cds.healpix.Projection.Y_INDEX;

import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.FOUR_OVER_PI;

import static cds.healpix.common.math.HackersDelight.toBits;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
import static cds.healpix.common.math.FastMath.sinQ;
import static cds.healpix.common.math.HackersDelight.BUT_SIGN_BIT_MASK_L;

import static cds.healpix.Healpix.checkLatitude;
import static cds.healpix.Healpix.sqrtOfThreeTimeOneMinusSinOf;
import static cds.healpix.Healpix.TRANSITION_Z;
import static cds.healpix.Healpix.TRANSITION_LATITUDE;
import static cds.healpix.Healpix.ONE_OVER_TRANSITION_Z;


/**
 * Convenience class made to manipulate temporary variables without having to create objects 
 * passed through method calls.
 * WARNING: the immediate consequence is that the class is not thread safe!
 * 
 * @author F.-X. Pineau
 *
 */
final class HealpixProjector {

  private double absLon;
  private long signLon;

  private double absLat;
  private long signLat;
  
  private int xOffset;
  
  private double x;
  private double y;
  
  
  /** Unable to instantiate out of the package. */
  HealpixProjector() { }
  
  /**
   * Make an HEALPix projection. The scale is such that base cell vertices coordinates are integers.
   * This projection is multi-purpose in the sense that:
   *  - if input lon in [-pi, pi], then output x in [-4, 4]
   *  - if input lon in [ 0, 2pi], then output x in [0, 8]
   *  - output y always in [-2, 2]
   * More generally
   *  - if input lon < 0, then output x in [-8, 0[
   *  - if input lon > 0, then output x in [ 0, 8[
   * WARNING: this method is not thread safe!!
   * @param lon longitude in radians, support positive and negative reasonably large values with a
   *        naive approach (no Cody-Waite nor Payne Hanek range reduction).
   * @param lat latitude, must be in [-pi/2, pi/2] radians
   * @param resultXY if lon <= 0, then x in [-8, 0], else if lon >= 0, x in [0, 8]
   *                 y always in [-2, 2].
   */
  void project(final double lon, final double lat, final double[] resultXY) {
    checkLatitude(lat);
    computeAbsOfLonLatAndKeepSigns(lon, lat);
    scaleAbsLonFromMinus1ToPlus1AndKeepRangeOffset();
    reduceRangeOffsetToMax7();
    if (isInEquatorialRegion()) {
      projectCylindricalEqualArea();
    } else {
      projectCollignon();
    }
    shiftXByReducedRangeOffset();
    applyLonLatSignsToXY();
    storeResultIn(resultXY);
  }

  private void computeAbsOfLonLatAndKeepSigns(final double lon, final double lat) {
    computeAbsOfLonAndKeepSign(lon);
    computeAbsOfLatAndKeepSign(lat);
  }
  
  private void computeAbsOfLonAndKeepSign(final double lon) {
    signLon = toBits(lon);                               assert lon == fromBits(signLon) : lon + " != " + fromBits(signLon);
    absLon = fromBits(signLon & BUT_SIGN_BIT_MASK_L);    assert  0 <= absLon && absLon <= 10 * TWO_PI;
    signLon &= SIGN_BIT_MASK_L;                          assert signLon == 0 || signLon == SIGN_BIT_MASK_L;
  }
  
  private void computeAbsOfLatAndKeepSign(final double lat) {
    signLat = toBits(lat);                               assert lat == fromBits(signLat);
    absLat = fromBits(signLat & BUT_SIGN_BIT_MASK_L);    assert  0 <= absLat && absLat <= HALF_PI; 
    signLat &= SIGN_BIT_MASK_L;                          assert signLat == 0 || signLat == SIGN_BIT_MASK_L;
  }

  private void scaleAbsLonFromMinus1ToPlus1AndKeepRangeOffset() {
    x = FOUR_OVER_PI * absLon;
    xOffset = ((int) x) | 1;      assert 1 <= xOffset && xOffset % 2 == 1;
    x -= xOffset;                 assert -1 <= x && x <= 1;
  }
  
  private void reduceRangeOffsetToMax7() { // <=> xOffset = xOffset modulo 8
    xOffset &= 7;    assert 1 <= xOffset && xOffset <= 7 && xOffset % 2 == 1;
  }
  
  private boolean isInEquatorialRegion() {
    return absLat <= TRANSITION_LATITUDE;
  }
  
  private void projectCylindricalEqualArea() {
    y = sinQ(absLat, 0);          assert  0 <= y && y <= TRANSITION_Z;
    y *= ONE_OVER_TRANSITION_Z;   assert  0 <= 0 && y <= 1;
  }
  
  private void projectCollignon() {
    y = sqrtOfThreeTimeOneMinusSinOf(absLat);    assert  0 <= y && y <  1; 
    x *= y;                                      assert -1 <= x && x <= 1;
    y = (2 - y);                                 assert  1 <  y && y <= 2;
  }

  private void shiftXByReducedRangeOffset() {
    x += xOffset;    assert  0 <= x && x <= 8;
  }
  
  private void applyLonLatSignsToXY() {
    applyLonSignToX();
    applyLatSignToY();
  }
  
  private void applyLonSignToX() {
    x = fromBits(signLon | toBits(x));    assert -8 < x && x < 8;
  }

  private void applyLatSignToY() {
    y = fromBits(signLat | toBits(y));    assert -2 <= y && y <= 2;
  }
 
  private void storeResultIn(final double[] resultXY) {
    resultXY[X_INDEX] = x;
    resultXY[Y_INDEX] = y;
  }

}
