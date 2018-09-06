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
import static cds.healpix.common.math.HackersDelight.floorIntP;
import static cds.healpix.common.math.HackersDelight.floorLongP;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.toBits;

import cds.healpix.fillingcurve.ZOrderCurve3D;

/**
 * Implementation of the HEALPix plus auxiliary axis indexation.
 * 
 * @author F.-X. Pineau
 *
 */
class HealpixNestedHashComputerWithAux implements HashComputerWithAux {
  
  /** 19, largest possible depth we can store on a signed positive long
  * (4 bits for base cells + 3 bits per depth + 3 remaining bits). */
  public static final int DEPTH_MAX = 16;
  
  private final HealpixNested h;
  private final AuxiliaryAxis auxAxis;
  private final long halfNside4IEEEdouble;
  
  private final HealpixProjector proj = new HealpixProjector();
  private double[] xy = new double[2];  // Result of the HEALPix projection
  private double x, y;                  // Coordinates in the shifted, rotated, scaled frame
  private long xInt, yInt;              // Cell coordinates in the shifted, rotated, scaled frame
  private int iBaseCell, jBaseCell;     // Base cell coos in the shifted, rotated, re-scaled frame
  private int iInBaseCell, jInBaseCell; // South-East / South-West coordinates inside a base cell
  private long baseCellBits;             // Hash value of the base cell (i.e. the depth 0 cell), shifted

  HealpixNestedHashComputerWithAux(final HealpixNested healpixNested, final AuxiliaryAxis auxAxis) {
    this.h = healpixNested;
    this.auxAxis = auxAxis;
    if (this.h.depth > 19) {
      throw new IllegalArgumentException("depth with auxiliary axis limited to 19 (inclusive)");
    }
    this.halfNside4IEEEdouble = h.depth == 0 ? 0L : (h.depth - 1L) << 52;
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
  public long hash(final double lonRad, final double latRad, final double auxValue) {
    project(lonRad, latRad);
    shiftAndRotateAndScale();
    discretize();
    computeBaseCellCoos();
    computeCoosInBaseCell();
    computeBaseCellBits();
    int aux = computeAuxValue(auxValue);
    return baseCellBits | ZOrderCurve3D.INSTANCE.ijk2hash(iInBaseCell, jInBaseCell, aux);
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
    x += y;
    y -= xy[X_INDEX];
    // Shift along the new y-axis so all the full projection contains only positive coordinates
    // to avoid a floor on both positive and negative values (to avoid a branch)
    y += 8;
    // Scale (multiply by nside/2) so that base-resolution cells are of size nside x nside
    x = fromBits(this.halfNside4IEEEdouble + toBits(x));
    y = fromBits(this.halfNside4IEEEdouble + toBits(y));
    assert 0 <= x && x <= 5 * h.nside;
    assert 0 <= y && y <= 5 * h.nside;
  }

  private void discretize() {
    xInt = floorLongP(x);    assert 0 <= xInt && xInt <= 5 * h.nside : xInt + " <= " + 4 * h.nside;
    yInt = floorLongP(y);    assert 0 <= yInt && yInt <= 5 * h.nside : yInt + " <= " + 4 * h.nside;
  }

  private void computeBaseCellCoos() {
    iBaseCell = h.dividedByNsideQuotient(xInt);    assert 0 <= iBaseCell && iBaseCell <= 4;
    jBaseCell = h.dividedByNsideQuotient(yInt);    assert 0 <= jBaseCell && jBaseCell <= 4;
  }

  private void computeCoosInBaseCell() {
    iInBaseCell = h.moduloNside(xInt);    assert 0 <= iInBaseCell && iInBaseCell < h.nside;
    jInBaseCell = h.moduloNside(yInt);    assert 0 <= jInBaseCell && jInBaseCell < h.nside;
  }
  
  private int computeAuxValue(final double auxValue) {
    final double normalized = this.auxAxis.project(auxValue);
    if (normalized < 0 || 1 <= normalized) {
      throw new IllegalArgumentException("Wrong auxiliary value: f(auxValue) must be in [0, 1[");
    }
    return floorIntP(normalized * this.h.nside);
    
  }
  
  /*private void computeBaseCellBits() {
    baseCellBits = h.d0cBitsLUPT[iBaseCell][jBaseCell] << h.depth;
  }*/
  
  private void computeBaseCellBits() {
    jBaseCell = 5 - (iBaseCell + jBaseCell);
    if (jBaseCell >= 0) {        
      assert jBaseCell <= 2;
      baseCellBits = ((long) ((jBaseCell << 2) + ((iBaseCell - ((--jBaseCell) >>> 63)) & 3))) << (h.twiceDepth + h.depth);
    } else if (jBaseCell == -1) { // rare, so few risks of branch miss-prediction
      baseCellBits = ((((long) ((iBaseCell - 1) & 3))) << (h.twiceDepth + h.depth)) | h.yMask;
    } else if (jBaseCell == -2) { // rare, so few risks of branch miss-prediction
      baseCellBits = (((long) (iBaseCell - 2)) << (h.twiceDepth + h.depth)) | h.xyMask;
    } else {
      assert false;
    }
  }
  
}
