package cds.healpix;

import static cds.healpix.Healpix.getBestStartingDepth;
import static cds.healpix.HealpixNestedBMOC.buildValue;
import static cds.healpix.NestedSmallCell.normalizeLon;
import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;

import java.util.Arrays;
import java.util.EnumSet;

import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.common.sphgeom.EllipticalCone;


public class NestedEllipticalConeComputerApprox {

  private static final EnumSet<Cardinal> ALL_CARDINALS = EnumSet.allOf(Cardinal.class);
  private static final double SQRT3 = Math.sqrt(3);

  private final int startDepth;
  private final HealpixNested startHpx; 
  private final HealpixNested deepHpx;
  private final VerticesAndPathComputer[] hcc;

  private final double aRad;
  private final double bRad;
  private final double paRad;
  private final double sinasinb;
  private EllipticalCone ellipse;
  
  private final NeighbourSelector neigSelector;
  private final FlatHashList neigList;

  private final int deltaDepthMax;

  public static enum Mode {
    OVERLAPPING_CELLS() {
      @Override
      public boolean isOk(final EllipticalCone ellipse,
          final VerticesAndPathComputer vpc, final long hash, 
          final double cellCenterLon, final double cellCenterLat) {
        return true;
      }
    },
    OVERLAPPING_CENTERS() {
      @Override
      public boolean isOk(final EllipticalCone ellipse,
          final VerticesAndPathComputer vpc, final long hash, 
          final double cellCenterLon, final double cellCenterLat) {
        return ellipse.contains(cellCenterLon, cellCenterLat);
      }
    },
    FULLY_IN() {
      @Override
      public boolean isOk(final EllipticalCone ellipse,
          final VerticesAndPathComputer vpc, final long hash, 
          final double cellCenterLon, final double cellCenterLat) {
        for (final double[] vertex : vpc.vertices(hash, ALL_CARDINALS).values()) {
          if (!ellipse.contains(vertex[LON_INDEX], vertex[LAT_INDEX])) {
            return false;
          }
        }
        return true;
      }
    };
    public abstract boolean isOk(EllipticalCone ellipse, VerticesAndPathComputer vpc,
        long hash, final double cellCenterLon, final double cellCenterLat);
  }

  /**
   * 
   * @param aRad elliptical cone major axis
   * @param hDeepestCells 
   */
  public NestedEllipticalConeComputerApprox(final double aRad, final double bRad, final double posAngRad,
      final HealpixNested deepHpx) {
    this.deepHpx = deepHpx;
    this.aRad = aRad;
    this.bRad = bRad;
    this.paRad = posAngRad;
    this.ellipse = new EllipticalCone(0, 0, aRad, bRad, posAngRad);
    this.sinasinb = this.ellipse.getSinA() * this.ellipse.getSinB();
    final int optimalStartingDepth = Math.min(this.deepHpx.depth, getBestStartingDepth(aRad));
    if (optimalStartingDepth == -1) {
      this.startDepth = 0;
      this.startHpx = null;
      this.neigSelector = null;
      this.deltaDepthMax = this.deepHpx.depth;
    } else {
      this.startDepth = optimalStartingDepth;
      this.startHpx = Healpix.getNested(this.startDepth);
      this.neigSelector = this.startHpx.newNeighbourSelector();
      this.deltaDepthMax = this.deepHpx.depth - this.startDepth;
    }
    this.neigList = new FlatHashList(-1, 9); // We don't care about the depth, internal usage only

    this.hcc = new VerticesAndPathComputer[this.deltaDepthMax + 1];
    for (int i = 0; i <= this.deltaDepthMax; i++) {
      this.hcc[i] = Healpix.getNested(this.startDepth + i).newVerticesAndPathComputer();
    }
  }

  public HealpixNestedBMOC overlapping(double ellipseCenterLonRad, double ellipseCenterLatRad,
      final Mode mode) {
    //this.ellipse.setProjCenter(ellipseCenterLonRad, ellipseCenterLatRad);
    this.ellipse = new EllipticalCone(ellipseCenterLonRad, ellipseCenterLatRad, aRad, bRad, paRad);
    // Store required space in the MOC
    final long[] mocElems = new long[this.nMocCellInAreaUpperBound()];
    int mocSize = 0;
    if (this.startHpx == null) {
      for (int h = 0; h < 12; h++) {
        mocSize = buildMocRecursively(mocElems, mocSize, 0, h, mode);
      }
    } else {
      // Compute hash of the cell containing the cone center
      final long centerHash = this.startHpx.hash(ellipseCenterLonRad, ellipseCenterLatRad); // lat is checked here!
      this.neigSelector.neighbours(centerHash, this.neigList);
      this.neigList.put(centerHash);
      this.neigList.sortByHashAsc();
      for (int i = 0; i < this.neigList.size(); i++) {
        mocSize = buildMocRecursively(mocElems, mocSize, 0, this.neigList.get(i), mode);
      }
    }
    return HealpixNestedBMOC.createPacking(this.deepHpx.depth, mocElems, mocSize);
  }

  private final int buildMocRecursively(final long[] moc, int mocLength, int deltaDepth, long hash,
      final Mode mode) {
    final int depth = this.startDepth + deltaDepth;
    // System.out.println("depth: " + depth);
    assert this.hcc[deltaDepth].depth() == depth : this.hcc[deltaDepth].depth() + " != " + depth;
    final VerticesAndPathComputer vpc = this.hcc[deltaDepth];
    final double[] center = vpc.center(hash);
    final double cellCenterLon = center[LON_INDEX];
    final double cellCenterLat = center[LAT_INDEX];
    final double rCircumCircle = Healpix.getLargestCenterToCellVertexDistance(
        cellCenterLon, cellCenterLat, depth);
    if (this.ellipse.containsCone(cellCenterLon, cellCenterLat, rCircumCircle)) {
      // we could have called ellipse.contains on the 4 cell vertices (more precise but time consuming)
      moc[mocLength++] = buildValue(depth, hash, true, this.deepHpx.depth);
      // System.out.println("add depth: " + depth + "; hash: " + hash + "; circumDeg: " + Math.toDegrees(rCircumCircle));
    } else if (this.ellipse.overlapCone(cellCenterLon, cellCenterLat, rCircumCircle)) {
      if (depth == this.deepHpx.depth) {
        if (mode.isOk(this.ellipse, vpc, hash, cellCenterLon, cellCenterLat)) { 
          moc[mocLength++] = buildValue(depth, hash, false, this.deepHpx.depth);
        }
      } else {
        hash <<= 2;
        deltaDepth++;
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth,   hash, mode);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, mode);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, mode);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, mode);
      }
    } // else cell fully out of the cone
    return mocLength;
  }

  private int nMocCellInAreaUpperBound() {
    // cell_area = 4 * pi / ncell = 4 * pi / (3 * 4 * nside^2) = pi / (3 * nside^2) =  pi * r^2
    // cell_radius = r = 1 / (sqrt(3) * nside)
    // As a very simple and naive rule, we take 4x the number of cells needed to cover
    // the cone external annulus
    // Annulus area = pi (sin(a+r)sinn(b+r) - sina*sinb)
    // N cells = 4 * pi (sin(a+r)sinn(b+r) - sina*sinb) /pi r^2 
    //         = 4 * 3 * nside^2 * (sin(a+r)*sin(b+r) - sina*sinb)
    final double oneOverR2 = 3 * (this.deepHpx.nside * this.deepHpx.nside);
    final double r = 1 / (SQRT3 * this.deepHpx.nside);
    return (int) (this.deepHpx.nHash * (1 + ((this.ellipse.getA()  + r) * (this.ellipse.getB()  + r) - this.ellipse.getA() * this.ellipse.getB())));
  }
}
