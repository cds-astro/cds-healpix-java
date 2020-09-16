package cds.healpix;

import static cds.healpix.HealpixNested.bits2hash;
import static cds.healpix.Projection.X_INDEX;
import static cds.healpix.Projection.Y_INDEX;
import static cds.healpix.common.math.HackersDelight.floorLongP;
import static cds.healpix.common.math.HackersDelight.toBits;

import java.util.EnumSet;

import cds.healpix.CompassPoint.MainWind;
import static cds.healpix.CompassPoint.MainWind.*;


/**
 *  WARNING: not tested yet!
 *  
 * @author F.-X. Pineau
 *
 */
class HealpixNestedInterpolator implements Interpolator {
  
  private final HealpixNested h;
  private final NeighbourSelector hns;

  private final long EQUAT_D0H_BITS_MASK;
  
  private final HealpixProjector proj = new HealpixProjector();
  private double[] xy = new double[2];  // Result of the HEALPix projection
  private double x, y;                  // Coordinates in the shifted, rotated, scaled frame
  private long xInt, yInt;              // Cell coordinates in the shifted, rotated, scaled frame
  private int iBaseCell, jBaseCell;     // Base cell coos in the shifted, rotated, re-scaled frame
  private int iInBaseCell, jInBaseCell; // South-East / South-West coordinates inside a base cell
  private long baseCellBits;            // Hash value of the base cell (i.e. the depth 0 cell), shifted

  private double cellCenterX, cellCenterY;
  private double dX, dY;
  private int quarter = 0; // 0 => S, 1 => E, 2 => W, 3 => N

  private NeighbourList neigList;
  
  HealpixNestedInterpolator(final HealpixNested healpixNested) {
    this.h = healpixNested;
    this.hns = new HealpixNestedNeighbourSelector(this.h); // thread safe
    this.neigList = new NeighbourList(this.h.depth);
    this.EQUAT_D0H_BITS_MASK = 4L << h.twiceDepth;
  }

  @Override
  public int depth() {
    return h.depth;
  }
  
  @Override
  public int interpolate(double lonRad, double latRad, long[] resultHashs, double[] resultWeights) {
    // Hash computation
    project(lonRad, latRad);
    shiftAndRotateAndScale();
    discretize();
    computeQuarter(); // in addition to simple hash computation
    computeBaseCellCoos();
    computeCoosInBaseCell();
    computeBaseCellHashBits();
    return computeHashsAndWeights(resultHashs, resultWeights);
  }
  
  /*private void computeWeigthsHashAndWeight(int dx, int dy, int i,
      long[] resultHashs, double[] resultWeights) {
    resultHashs[i] = baseCellBits | h.fc.ij2hash(iInBaseCell + dx, jInBaseCell + dy);
    resultWeights[i] = (dx + dX) * (dy + dY);
  }*/
  
  private boolean isBaseCellInEquatRegionFromBits(final long baseCellBits) {
    return (baseCellBits & EQUAT_D0H_BITS_MASK) != 0L;
  }
  
  private boolean isInNorthPolarCapFromBits(final long baseCellBits) {
    return baseCellBits > EQUAT_D0H_BITS_MASK;
  }
  
  private int computeHashsAndWeights(long[] resultHashs, double[] resultWeights) {
    final long ijInterleavedBits = h.fc.ij2hash(iInBaseCell, jInBaseCell);
    final long hashC = baseCellBits | ijInterleavedBits;
    
    if (isInBaseCellBorderFromBits(iInBaseCell, jInBaseCell)) { // border, many case to handle :o/
      // We accept to be slower for those 24 cells in the sphere
      // - we introduce redundancy, not to have to add the "neighbours" code here!
      if (isBaseCellInEquatRegionFromBits(baseCellBits)) {
        // if S or N => no N or S neigh
        if (ijInterleavedBits == 0L && quarter == 0) { // South 
          // Add a fake cell having for sum the mean of SE + SW
          // => C * w_c + SE * w_se + SW * w_sw + (SE + SW) / 2 * w_s
          // => C * w_c + SE * (w_se + 1/2 * w_s) + SW * (w_sw + 1/2 * w_s)
          assert dX >= 0 && dY >= 0;
          // C
          resultHashs[0] = hashC;
          resultWeights[0] = (1 - dX) * (1 - dY);
          double halfWS = 0.5 * dX * dY;
          // SW
          resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell - 1, jInBaseCell);
          resultWeights[2] = dX * (1 - dY) + halfWS;
          // SE
          resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell - 1);
          resultWeights[3] = (1 - dX) * dY + halfWS;
          return 3;
        } else if (ijInterleavedBits == h.xyMask && quarter == 3) { // North
          assert dX < 0 && dY < 0;
          // C
          resultHashs[0] = hashC;
          resultWeights[0] = (1 + dX) * (1 + dY);
          // N
          double halfWN = 0.5 * dX * dY;
          // NE
          resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell + 1, jInBaseCell);
          resultWeights[2] = dX * (1 + dY) + halfWN;
          // NW
          resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell + 1);
          resultWeights[3] = (1 + dX) * dY + halfWN;
          return 3;
        } else { // Neither S nor N corner
          return computeHashsAndWeights(hashC, quarter, resultHashs, resultWeights); 
        }
      } else if (isInNorthPolarCapFromBits(baseCellBits)) { // north polar cap
        // if E or W => no E or W neigh
        if (ijInterleavedBits == h.xMask && quarter == 1) { // East
          assert dX < 0 && dY >= 0;
          // E
          double halfWE = -0.5 * dX * dY; // to be added to C and NE
          // C
          resultHashs[0] = hashC;
          resultWeights[0] = (1 + dX) * (1 - dY) + halfWE;
          // NE
          resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell + 1, jInBaseCell);
          resultWeights[2] = -dX * (1 - dY) + halfWE;
          // SE
          resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell - 1);
          resultWeights[3] = (1 - dX) * dY;
          return 3;
        } else if (ijInterleavedBits == h.yMask && quarter == 2) { // West
          assert dX >= 0 && dY < 0;
          // W
          double halfWW = -0.5 * dX * dY; // to be added to C and NW
          // C
          resultHashs[0] = hashC;
          resultWeights[0] = (1 - dX) * (1 + dY) + halfWW;
          // SW
          resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell - 1, jInBaseCell);
          resultWeights[2] = dX * (1 + dY);
          // NW
          resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell + 1);
          resultWeights[3] = (1 - dX) * dY + halfWW;
          return 3;
        } else { // Neither E nor W corner
          return computeHashsAndWeights(hashC, quarter, resultHashs, resultWeights); 
        }
      } else { // south polr cap
        // if E or W => no E or W neigh
        if (ijInterleavedBits == h.xMask && quarter == 1) { // East
          assert dX < 0 && dY >= 0;
          // E
          double halfWE = -0.5 * dX * dY; // to be added to C and SE
          // C
          resultHashs[0] = hashC;
          resultWeights[0] = (1 + dX) * (1 - dY) + halfWE;
          // NE
          resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell + 1, jInBaseCell);
          resultWeights[2] = -dX * (1 - dY);
          // SE
          resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell - 1);
          resultWeights[3] = (1 - dX) * dY + halfWE;
          return 3;
        } else if (ijInterleavedBits == h.yMask && quarter == 2) { // West
          assert dX >= 0 && dY < 0;
          // W
          double halfWW = -0.5 * dX * dY; // to be added to C and SW
          // C
          resultHashs[0] = hashC;
          resultWeights[0] = (1 - dX) * (1 + dY) + halfWW;
          // SW
          resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell - 1, jInBaseCell);
          resultWeights[2] = dX * (1 + dY) + halfWW;
          // NW
          resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell + 1);
          resultWeights[3] = (1 - dX) * dY;
          return 3;
        } else { // Neither E nor W corner
          return computeHashsAndWeights(hashC, quarter, resultHashs, resultWeights); 
        }
      }
    } else { // easy, all cells except 24
      // It would have been cleaner to call here the extern neighbours method.
      // We introduce here redundancy for perf purposes (at the cost of duplication on redability)
      switch(quarter) {
      case 0: // S
        assert dX >= 0 && dY >= 0;
        // C
        resultHashs[0] = hashC;
        resultWeights[0] = (1 - dX) * (1 - dY);
        // S
        resultHashs[1] = baseCellBits | h.fc.ij2hash(iInBaseCell - 1, jInBaseCell - 1);
        resultWeights[1] = dX * dY;
        // SW
        resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell - 1, jInBaseCell);
        resultWeights[2] = dX * (1 - dY);
        // SE
        resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell - 1);
        resultWeights[3] = (1 - dX) * dY;
        break;
      case 1: // E
        assert dX < 0 && dY >= 0;
        // C
        resultHashs[0] = hashC;
        resultWeights[0] = (1 + dX) * (1 - dY);
        // E
        resultHashs[1] = baseCellBits | h.fc.ij2hash(iInBaseCell + 1, jInBaseCell - 1);
        resultWeights[1] = -dX * dY;
        // NE
        resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell + 1, jInBaseCell);
        resultWeights[2] = -dX * (1 - dY);
        // SE
        resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell - 1);
        resultWeights[3] = (1 - dX) * dY;
        break;
      case 2: // W
        assert dX >= 0 && dY < 0;
        // C
        resultHashs[0] = hashC;
        resultWeights[0] = (1 - dX) * (1 + dY);
        // W
        resultHashs[1] = baseCellBits | h.fc.ij2hash(iInBaseCell - 1, jInBaseCell + 1);
        resultWeights[1] = -dX * dY;
        // SW
        resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell - 1, jInBaseCell);
        resultWeights[2] = dX * (1 + dY);
        // NW
        resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell + 1);
        resultWeights[3] = (1 - dX) * dY;
        break;
      case 3: // N
        assert dX < 0 && dY < 0;
        // C
        resultHashs[0] = hashC;
        resultWeights[0] = (1 + dX) * (1 + dY);
        // N
        resultHashs[1] = baseCellBits | h.fc.ij2hash(iInBaseCell + 1, jInBaseCell + 1);
        resultWeights[1] = dX * dY;
        // NE
        resultHashs[2] = baseCellBits | h.fc.ij2hash(iInBaseCell + 1, jInBaseCell);
        resultWeights[2] = dX * (1 + dY);
        // NW
        resultHashs[3] = baseCellBits | h.fc.ij2hash(iInBaseCell, jInBaseCell + 1);
        resultWeights[3] = (1 + dX) * dY;
        break;
      default:
        assert false;
      }
      return 4;
    }
  }
  
  private int computeHashsAndWeights(final long hashC, int quarter, 
      long[] resultHashs, double[] resultWeights) {
    switch(quarter) {
    case 0: // S
      assert dX >= 0 && dY >= 0;
      this.hns.neighbours(hashC, EnumSet.of(S, SE, SW), this.neigList);
      // C
      resultHashs[0] = hashC;
      resultWeights[0] = (1 - dX) * (1 - dY);
      // S
      resultHashs[1] = this.neigList.get(S);
      resultWeights[1] = dX * dY;
      // SW
      resultHashs[2] = this.neigList.get(SW);
      resultWeights[2] = dX * (1 - dY);
      // SE
      resultHashs[3] = this.neigList.get(SE);
      resultWeights[3] = (1 - dX) * dY;
      break;
    case 1: // E
      assert dX < 0 && dY >= 0;
      this.hns.neighbours(hashC, EnumSet.of(E, SE, SW), this.neigList);
      // C
      resultHashs[0] = hashC;
      resultWeights[0] = (1 + dX) * (1 - dY);
      // E
      resultHashs[1] = this.neigList.get(E);
      resultWeights[1] = -dX * dY;
      // NE
      resultHashs[2] = this.neigList.get(NE);
      resultWeights[2] = -dX * (1 - dY);
      // SE
      resultHashs[3] = this.neigList.get(SE);
      resultWeights[3] = (1 - dX) * dY;
      break;
    case 2: // W
      assert dX >= 0 && dY < 0;
      this.hns.neighbours(hashC, EnumSet.of(W, SW, NW), this.neigList);
      // C
      resultHashs[0] = hashC;
      resultWeights[0] = (1 - dX) * (1 + dY);
      // W
      resultHashs[1] = this.neigList.get(W);
      resultWeights[1] = -dX * dY;
      // SW
      resultHashs[2] = this.neigList.get(SW);
      resultWeights[2] = dX * (1 + dY);
      // NW
      resultHashs[3] = this.neigList.get(NW);
      resultWeights[3] = (1 - dX) * dY;
      break;
    case 3: // N
      assert dX < 0 && dY < 0;
      this.hns.neighbours(hashC, EnumSet.of(N, NE, NW), this.neigList);
      // C
      resultHashs[0] = hashC;
      resultWeights[0] = (1 + dX) * (1 + dY);
      // N
      resultHashs[1] = this.neigList.get(N);
      resultWeights[1] = dX * dY;
      // NE
      resultHashs[2] = this.neigList.get(NE);
      resultWeights[2] = dX * (1 + dY);
      // NW
      resultHashs[3] = this.neigList.get(NW);
      resultWeights[3] = (1 + dX) * dY;
      break;
    default:
      assert false;
    }
    this.neigList.clear();
    return 4;
  }
  
  private boolean isInBaseCellBorderFromBits(
      final long iInBaseCell, final long jInBaseCell) {
    return 0 == iInBaseCell || iInBaseCell == this.h.nsideRemainderMask
        || 0 == jInBaseCell || jInBaseCell == this.h.nsideRemainderMask;
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
    x = h.timeHalfNside(x);
    y = h.timeHalfNside(y);
    assert 0 <= x && x <= (5+1e-15) * h.nside : x;
    assert 0 <= y && y <= (5+1e-15) * h.nside : y;
  }

  private void discretize() {
    xInt = floorLongP(x);    assert 0 <= xInt && xInt <= 5 * h.nside : xInt + " <= " + 4 * h.nside;
    yInt = floorLongP(y);    assert 0 <= yInt && yInt <= 5 * h.nside : yInt + " <= " + 4 * h.nside;
  }

  private void computeQuarter() {
    // Watning: special treatment to be done if border!!
    cellCenterX = xInt + 0.5d;
    cellCenterY = yInt + 0.5d;
    
    dX = cellCenterX - x;
    dY = cellCenterY - y;
    
    quarter  = (int) (toBits(dX) >> 63); // +0 OR +1
    assert (x >= cellCenterX && quarter == 1) || (x < cellCenterX && quarter == 0);
    
    quarter |= (int) (toBits(dY) >> 62); // +0 or +2
    assert (y >= cellCenterY && (quarter | 2) == 2) || (y < cellCenterY && (quarter | 2) == 0);

    assert quarter >= 0 && quarter < 4;
  }
  
  private void computeBaseCellCoos() {
    iBaseCell = h.dividedByNsideQuotient(xInt);    assert 0 <= iBaseCell && iBaseCell <= 5;
    jBaseCell = h.dividedByNsideQuotient(yInt);    assert 0 <= jBaseCell && jBaseCell <= 5 : jBaseCell;
  }

  private void computeCoosInBaseCell() {
    iInBaseCell = h.moduloNside(xInt);    assert 0 <= iInBaseCell && iInBaseCell < h.nside;
    jInBaseCell = h.moduloNside(yInt);    assert 0 <= jInBaseCell && jInBaseCell < h.nside;
  }
  /**
   Depending on the hardware, and the soft cache occupation, better or not than the other version
  */
  /*private void computeBaseCellHashBits() {
    baseCellBits = h.d0cBitsLUPT[iBaseCell][jBaseCell];    assert baseCellBits >= 0;
    assert Math.abs((4 - jBaseCell) - iBaseCell) < 2
      || (iInBaseCell == 0 && jInBaseCell == 0 && ((baseCellBits & h.xyMask) == h.xyMask))
      || (jInBaseCell == 0 && ((baseCellBits & h.yMask) == h.xyMask));
  }*/
  
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
  }

}
