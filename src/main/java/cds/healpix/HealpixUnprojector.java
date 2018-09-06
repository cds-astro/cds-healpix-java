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
import static cds.healpix.Healpix.TRANSITION_LATITUDE;

import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;
import static cds.healpix.common.math.HackersDelight.toBits;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
import static cds.healpix.common.math.FastMath.acos;
import static cds.healpix.common.math.FastMath.asin;
import static cds.healpix.common.math.HackersDelight.BUT_SIGN_BIT_MASK_L;

import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;

final class HealpixUnprojector {

  /** 1/sqrt(6), constant used in the Collignon deprojection. */
  public static final double ONE_OVER_SQRT6 = 0.40824829046386301636;
  
  static final double EPS_POLE = 1e-13;
  static final double EPS_EDGE = 1e-12;
  
  
  private long signX;
  private double absX;
  
  private long signY;
  private double absY;
  
  private int xOffset;
  
  private double lon;
  private double lat;
  
  /** Unable to instantiate out of the package. */
  HealpixUnprojector() { }
  
  /**
   * Unproject HEALPix projected points.
   * This unprojection is multi-purpose in the sense that:
   *  - if input x in [-8, 0[, then output lon in [-2pi, 0]
   *  - if input x in [ 0, 8], then output x in [0, 2pi]
   *  - output lat always in [-pi/2, pi/2]
   * WARNING: this method is not thread safe!!
   * @param x x coordinate, support positive and negative reasonably large values with a
   *        naive approach (no Cody-Waite nor Payne Hanek range reduction).
   * @param y y coordinate, must be in [-2, 2]
   * @param resultLonLat if x <= 0, then lon in [-2pi, 0], else if x >= 0, lon in [0, 2pi]
   *                     lat always in [-pi/2, pi/2].
   */
  void unproject(final double x, final double y, final double[] resultLonLat) {
    checkProjectionBoundsInY(y);
    computeAbsOfXYAndKeepSigns(x, y);
    shiftAbsXInMinus1ToPlus1RangeAndKeepRangeOffset();
    reduceRangeOffsetToMax7();
    if (isInEquatorialRegion()) {
      deprojCylindricalEquaArea(y);
    } else {
      deprojCollignon(x);
    }
    shiftLonToReducedRangeAndScaleToRadians();
    applyXYSignsToLonLat();
    storeResultIn(resultLonLat);
  }
  
  private static void checkProjectionBoundsInY(final double y) {
    if (y < -2 || y > 2) {
      throw new IllegalArgumentException("Wrong y. Expected: in [-2, 2]. Actual: " + y);
    }
  }
  
  private void computeAbsOfXYAndKeepSigns(final double x, final double y) {
    computeAbsOfXAndKeepSign(x);
    computeAbsOfYAndKeepSign(y);
  }
  
  private void computeAbsOfXAndKeepSign(final double x) {
    signX = toBits(x);                               assert x == fromBits(signX);
    absX = fromBits(signX & BUT_SIGN_BIT_MASK_L);    assert 0 <= absX;
    signX &= SIGN_BIT_MASK_L;                        assert 0 == signX || SIGN_BIT_MASK_L == signX;
  }
  
  private void computeAbsOfYAndKeepSign(final double y) {
    signY = toBits(y);                               assert y == fromBits(signY);
    absY = fromBits(signY & BUT_SIGN_BIT_MASK_L);    assert 0 <= absY && absY <= 2; 
    signY &= SIGN_BIT_MASK_L;                        assert 0 == signY || SIGN_BIT_MASK_L == signY;
  }
  
  private void shiftAbsXInMinus1ToPlus1RangeAndKeepRangeOffset() {
    xOffset = ((int) absX) | 1;    assert  1 <= xOffset && xOffset % 2 == 1;
    lon = absX - xOffset;          assert -1 <= lon && lon <= 1;
  }
  
  private void reduceRangeOffsetToMax7() {
    xOffset &= 7;    assert  1 <= xOffset && xOffset <= 7 && xOffset % 2 == 1;
  }
  
  private boolean isInEquatorialRegion() {
    return absY <= 1;
  }
  
  private void deprojCylindricalEquaArea(final double y) {
    lat = asin(y * TRANSITION_Z);    assert -HALF_PI <= lat && lat <= HALF_PI;
  }
  
  private void deprojCollignon(final double x) {
    lat = 2 - absY;                       assert   0 <= lat && lat <  1;
    if (isNotNearFromPole(lat)) {      // Rare, so few risks of branch miss-prediction
      lon /= lat;                         assert  -1 - EPS_EDGE <= lon && lon <= 1 + EPS_EDGE;
      dealWithNumericalApproxInEdge();    assert  -1 <= lon && lon <= 1;                                 
    } // in case of pole, lon = lat = 0 (we avoid NaN due to division by lat=0)
    lat *= ONE_OVER_SQRT6;                assert   0 <= lat && lat < ONE_OVER_SQRT6;
    lat = 2 * acos(lat) - HALF_PI;        assert TRANSITION_LATITUDE < lat && lat <= HALF_PI;
  }
  
  private void dealWithNumericalApproxInEdge() { 
    // Both conditions are rare, so few risks of branch miss-prediction
    if (lon > 1) {            assert lon - 1 < EPS_EDGE;
      lon = 1;               
    } else if (lon < -1) {    assert lon + 1 > -EPS_EDGE;
      lon = -1;
    }
  } 
  
  static boolean isNotNearFromPole(final double sqrtOfThreeTimeOneMinusSinOf) {
    // In case of pole: x = y = 0
    return sqrtOfThreeTimeOneMinusSinOf > EPS_POLE;
  }
  
  private void shiftLonToReducedRangeAndScaleToRadians() {
    lon += xOffset;                       assert   0 <= lon && lon <= 8;
    lon *= PI_OVER_FOUR;                  assert   0 <= lon && lon <= TWO_PI;
  }
  
  private void applyXYSignsToLonLat() {
    applyXSignToLon();
    applyYSignToLat();
  }
  
  private void applyXSignToLon() {
    lon = fromBits(signX | toBits(lon));    assert -PI <= lat && lat <= TWO_PI;
  }
  
  private void applyYSignToLat() {
    lat = fromBits(signY | toBits(lat));    assert -HALF_PI <= lat && lat <= HALF_PI;
  }
  
  private void storeResultIn(final double[] resultLonLat) {
    resultLonLat[LON_INDEX] = lon;
    resultLonLat[LAT_INDEX] = lat;
  }

}
