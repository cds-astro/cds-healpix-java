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

import static cds.healpix.Healpix.checkLatitude;
import static cds.healpix.Healpix.sqrtOfThreeTimeOneMinusSinOf;
import static cds.healpix.Healpix.TRANSITION_LATITUDE;
import static cds.healpix.Healpix.ONE_OVER_TRANSITION_Z;
import static cds.healpix.Projection.X_INDEX;
import static cds.healpix.Projection.Y_INDEX;
import static cds.healpix.common.math.Math.FOUR_OVER_PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.SQRT6;

import cds.healpix.common.math.FastMath;

import static cds.healpix.common.math.HackersDelight.BUT_SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.floorLongP;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.toBits;

/**
 * This class simply prevents the creation of intermediary objects.
 * It must be used if you plan to call many times the {@link HealpixNested#hash(LonLat)} method.
 * WARNING: this class is not thread safe. Ask for one object of the type by thread in a
 * multithreaded environment.
 *
 * @author F.-X. Pineau
 *
 */
final class HealpixNestedHashComputer implements HashComputer {

  private final HealpixNested h;
  
  private double xpm1; // in [-1, 1]
  private int q;       // in [0, 3]: 0 if 0<=lon<pi/2; 1 if pi/2<=lon<pi/3; ...

  private long d0h;      // base cell (depth 0) hash value
  private double lInD0h; // x coordinate in the cartesian base cell diamon 
  private double hInD0h; // y coordiante in the cartesian base cell diamon
  
  HealpixNestedHashComputer(final HealpixNested healpixNested) {
    this.h = healpixNested;
  }

  @Override
  public int depth() {
    return h.depth;
  }

  /**
   * See {@link HashComputer#hash(double, double)}.
   * WARNING: in this implementation, the method is NOT THREAD-SAFE!
   */
  @Override
  public long hash(final double lonRad, final double latRad) {
    checkLatitude(latRad);
    this.d0hAndlhInD0c(lonRad, latRad);
    // Coords inside the base cell
    //  - ok to cast to int since small negative values due to numerical inaccuracies are rounded to 0    
    int i = (int) h.timeHalfNside(this.hInD0h + this.lInD0h);
    int j = (int) h.timeHalfNside(this.hInD0h - this.lInD0h);
    //  - deals with numerical inaccuracies, rare so branch miss-prediction negligible
    if (i == h.nside) { --i; } 
    if (j == h.nside) { --j; }
    // Build hash from components
    return (this.d0h << this.h.twiceDepth) | this.h.fc.ij2hash(i, j);
  }
  
  // Compute the base cell index (d0h) plus the h and l coordinates in it
  private void d0hAndlhInD0c(final double lonRad, final double latRad) {
    /// Transform the input longitude in a value `x` in `[-1, 1[` plus a quarter in `[0, 3]`,
    /// such that `lon = (x + 1) * PI / 4 + q * PI / 2`.
    this.xpm1AndQ(lonRad);
    assert -1 <= this.xpm1 && this.xpm1 <= 1;
    assert 0 <= this.q && this.q <= 3;
    /// Different branches for: north polar cap; south polar cap; equatorial region 
    if (latRad > TRANSITION_LATITUDE) {
      // North polar cap, Collignon projection.
      // - set the origin to (PI/4, 0)
      final double sqrt3OneMinZ = SQRT6 * FastMath.cosQ(0.5 * latRad + PI_OVER_FOUR);
      this.lInD0h = this.xpm1 * sqrt3OneMinZ;
      this.hInD0h = 2.0 - sqrt3OneMinZ;
      this.d0h = q;
    } else if (latRad < -TRANSITION_LATITUDE) {
      // South polar cap, Collignon projection
      // - set the origin to (PI/4, -PI/2)
      final double sqrt3OneMinZ = SQRT6 * FastMath.cosQ(PI_OVER_FOUR - 0.5 * latRad); // cos(-x) = cos(x)
      this.lInD0h = this.xpm1 * sqrt3OneMinZ;
      this.hInD0h = sqrt3OneMinZ;
      this.d0h = q + 8;
    } else {
      // Equatorial region, Cylindrical equal area projection
      // - set the origin to (PI/4, 0)               if q = 2
      // - set the origin to (PI/4, -PI/2)           if q = 0
      // - set the origin to (0, -TRANSITION_LAT)    if q = 3
      // - set the origin to (PI/2, -TRANSITION_LAT) if q = 1
      // All those lines because I took from FastMath sinQ, and sinQ take a positive angle
      long signLat = toBits(latRad);
      final double absLat = fromBits(signLat & BUT_SIGN_BIT_MASK_L);
      signLat &= SIGN_BIT_MASK_L;
      final double ypm1 = fromBits(toBits(FastMath.sinQ(absLat, 0) * ONE_OVER_TRANSITION_Z) | signLat);
      // Inequalities have been carefully chosen so that S->E and S->W axis are part of the cell,
      // and not E->N and W->N
      // Branch free version
      // |\2/|
      // .3X1.
      // |/0\|
      final int q01 = this.xpm1 >   ypm1 ? 1 : 0; /* 0/1 */ 
      final int q12 = this.xpm1 >= -ypm1 ? 1 : 0; /* 0\1 */ 
      final int q1 = q01 & q12; /* = 1 if q1, 0 else */
      final int q013 = q01 + (1 - q12); // = q01 + q03; /* 0/1 + 1\0 +  */
      // x: x_pm1 + 1 if q3 | x_pm1 - 1 if q1 | x_pm1 if q0 or q2
      this.lInD0h = this.xpm1 - ((q01 + q12) - 1);
      // y: y_pm1 + 0 if q2 | y_pm1 + 1 if q1 or q3 | y_pm1 + 2 if q0 
      this.hInD0h = ypm1 + q013;
      // d0h: +8 if q0 | +4 if q3 | +5 if q1
      this.d0h = (q013 << 2) + ((this.q + q1) & 3);
    }
  }
  
  /// Transform the input longitude, in radians, in a value `x` in `[-1, 1[` plus a quarter in `[0, 3]`,
  /// such that `lon = (x + 1) * PI / 4 + q * PI / 2`. 
  private void xpm1AndQ(double lonRad) { //(f64, u8) {
    // Compute absolute value of lon and keep its sign
    long lonSign = toBits(lonRad);
    final double lonAbs = fromBits(lonSign & BUT_SIGN_BIT_MASK_L);
    lonSign &= SIGN_BIT_MASK_L;
    // Compute x (scaled longitude so that [0, PI/2] -> [0, 2] and q, a number of PI/4 intervals
    final double x = lonAbs * FOUR_OVER_PI;
    final int q = ((int) x) | 1;  
    // Remark: to avoid the branch, we could have copied lon_sign on x - q, 
    //         but I so far lack of idea to deal with q efficiently.
    //         And we are not supposed to have negative longitudes in ICRS 
    //         (the most used reference system in astronomy).
    if (lonSign == 0L) { // => lon >= 0
      this.xpm1 = x - q;
      this.q = (q & 7) >> 1;
    } else { // case lon < 0, should be rare => few risks of branch miss-prediction
      // Since q is in [0, 3]: 3 - (q >> 1)) <=> 3 & !(q >> 1)
      // WARNING: BE SURE TO HANDLE THIS CORRECTLY IN THE REMAINING OF THE CODE!
      //  - Case lon =  3/4 pi = 270 deg => x = -1, q=3
      //  - Case lon = -1/2 pi = -90 deg => x =  1, q=2
      this.xpm1 = q - x;
      this.q = 3 - ((q & 7) >> 1);
    }
  }
  
  // The following code was quite elegant but suffer from the global projection numerical
  // inacurracies close to ra=n*pi/2 in the polar caps.
  // I though I dealed with those but a few corner cases remains.
  // Like in the Rust library, I prefer to replace this code with a new version I originally 
  // developped for graphic cards in Aladin Lite WebGL (then ported in glsl by Matthieu).
  // I (am supposed to have) improved the numerical precision working on (x, y, z) float inputs
  // and then re-adapted the code for (lon, lat) double inputs.
  // The corner cases are not present in HealpixNestedFast since we do not use the global 
  // projection but deal with polars caps and the eaquatorial region independently.
  
  /*private final HealpixProjector proj = new HealpixProjector();
  private double[] xy = new double[2];  // Result of the HEALPix projection
  private double x, y;                  // Coordinates in the shifted, rotated, scaled frame
  private long xInt, yInt;              // Cell coordinates in the shifted, rotated, scaled frame
  private int iBaseCell, jBaseCell;     // Base cell coos in the shifted, rotated, re-scaled frame
  private int iInBaseCell, jInBaseCell; // South-East / South-West coordinates inside a base cell
  private long baseCellBits;            // Hash value of the base cell (i.e. the depth 0 cell), shifted
  
  **
   * See {@link HashComputer#hash(double, double)}.
   * WARNING: in this implementation, the method is NOT THREAD-SAFE!
   *
  @Override
  public long hash(final double lonRad, final double latRad) {
    project(lonRad, latRad);
    shiftAndRotateAndScale();
    discretize();
    computeBaseCellCoos();
    computeCoosInBaseCell();
    computeBaseCellHashBits();
    return baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell);
  }

  private void project(final double lon, final double lat) {
    proj.project(lon, lat, xy);
    x = xy[X_INDEX];    assert -8 <  x && x <  8;
    y = xy[Y_INDEX];    assert -2 <= y && y <= 2;
    if (x < 0) { // I expect in most cases lon to be > 0, hence few chance of branch miss-prediction
      x += 8;    // else, could have written: x += (toBits(8) & SIGN_BIT_MASK_L) >> 60
    }
    assert 0 <= x && x <  8;
  }

  private void shiftAndRotateAndScale() {
    // Shift frame center from base cell 4 center to base cell 4 south vertex
    y++;
    // Rotate frame of -45 deg (cos and sin = +-\sqrt(2)/2) and scale by a factor sqrt(2)
    // so that the base-resolution cells are of size 2x2 with horizontal and vert. edges
    // x' = sqrt(2)/2 * ( cos * x + sin * y) * sqrt(2)
    // y' = sqrt(2)/2 * (-sin * x + cos * y) * sqrt(2)
    final double tmp = x;
    x += y;
    y -= tmp;
    // Shift along the new y-axis so all the full projection contains only positive coordinates
    // to avoid a floor on both positive and negative values (to avoid a branch)
    y += 8;
    // Scale (multiply by nside/2) so that base-resolution cells are of size nside x nside
    x = h.timeHalfNsideP(x);
    y = h.timeHalfNsideP(y);
    assert 0 <= x && x <= (5+1e-15) * h.nside : x;
    assert 0 <= y && y <= (5+1e-15) * h.nside : y;
  }

  private void discretize() {
    xInt = floorLongP(x);    assert 0 <= xInt && xInt <= 5L * h.nside : xInt + " <= " + 4L * h.nside;
    yInt = floorLongP(y);    assert 0 <= yInt && yInt <= 5L * h.nside : yInt + " <= " + 4L * h.nside;
  }

  private void computeBaseCellCoos() {
    iBaseCell = h.dividedByNsideQuotient(xInt);    assert 0 <= iBaseCell && iBaseCell <= 5;
    jBaseCell = h.dividedByNsideQuotient(yInt);    assert 0 <= jBaseCell && jBaseCell <= 5 : jBaseCell;
  }

  private void computeCoosInBaseCell() {
    iInBaseCell = h.moduloNside(xInt);    assert 0 <= iInBaseCell && iInBaseCell < h.nside;
    jInBaseCell = h.moduloNside(yInt);    assert 0 <= jInBaseCell && jInBaseCell < h.nside;
  }

  
  private void computeBaseCellHashBits() {
    jBaseCell = 5 - (iBaseCell + jBaseCell);
    if (jBaseCell >= 0) {        
      assert jBaseCell <= 2;
      baseCellBits = ((long) ((jBaseCell << 2) + ((iBaseCell - ((--jBaseCell) >>> 63)) & 3))) << h.twiceDepth;
    } else if (jBaseCell == -1) { // rare, so few risks of branch miss-prediction
      baseCellBits = ((((long) ((iBaseCell - 1) & 3))) << h.twiceDepth) | h.yMask;
    } else if (jBaseCell == -2) { // rare, so few risks of branch miss-prediction
      baseCellBits = (((long) (iBaseCell - 2)) << h.twiceDepth) | h.xyMask;
    } else {                      // should never enter this branch
      assert false;
    }
  }*/
  
}
