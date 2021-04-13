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

import static cds.healpix.HealpixNestedBMOC.buildValue;
import static cds.healpix.NestedSmallCell.normalizeLon;
import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.cos;

import java.util.logging.Logger;
import java.util.Arrays;
import java.util.EnumSet;

import cds.healpix.CompassPoint.Cardinal;

/**
 * WARNING: the flags in the BMOC returned by the method 'overlappingCells(double, double)'
 * may be aprroximation: a cell may be fully inside the cone, while its flag is set to 'PARTIAL'.
 *  
 * If you really want the list of cell fully inside, use the (slower) method
 * 'overlappingCells(double, double, FULLY_IN)'  
 * 
 * @author F.-X. Pineau
 *
 */
final class NestedSmallCellApproxedMethod implements HealpixNestedFixedRadiusConeComputer {

  private static final EnumSet<Cardinal> ALL_CARDINALS = EnumSet.allOf(Cardinal.class);

  private final AngularDistanceComputer angDistComputer;

  private final int startingDepth;
  private final int deeperDepth;
  private final int deltaDepthMax;

  private final double rRad;

  private final HashComputer hComputerStartingDepth;
  private final VerticesAndPathComputer[] hcc;
  private final NeighbourSelector neigSelector;

  private final HealpixNested hnDeeperDepth;
  private final HashComputer hComputerDeeperDepth;


  private final FlatHashList neigList;

  private int baseCellHash;             // Hash value of the base cell (i.e. the depth 0 cell)
  private int iInBaseCell, jInBaseCell; // South-East / South-West coordinates inside a base cell

  private static enum Mode {
    OVERLAPPING_CELLS() {
      @Override
      public boolean isOk(final double dConeCell, final double rRad,
          final VerticesAndPathComputer vpc, final long hash, 
          final double coneCenterLon, final double coneCenterLat, final double cosCenterLat,
          final AngularDistanceComputer angDistComputer) {
        return true;
      }
    },
    OVERLAPPING_CENTERS() {
      @Override
      public boolean isOk(final double dConeCell, final double rRad,
          final VerticesAndPathComputer vpc, final long hash, 
          final double coneCenterLon, final double coneCenterLat, final double cosCenterLat,
          final AngularDistanceComputer angDistComputer) {
        return dConeCell <= rRad;
      }
    },
    FULLY_IN() {
      @Override
      public boolean isOk(final double dConeCell, final double rRad,
          final VerticesAndPathComputer vpc, final long hash, 
          final double coneCenterLon, final double coneCenterLat, final double cosCenterLat,
          final AngularDistanceComputer angDistComputer) {
        return (dConeCell <= rRad) 
            && allVerticesOk(rRad, vpc, hash, coneCenterLon, coneCenterLat, cosCenterLat, angDistComputer);
      }
      private boolean allVerticesOk(final double rRad,
          final VerticesAndPathComputer vpc, final long hash, 
          final double coneCenterLon, final double coneCenterLat, final double cosCenterLat,
          final AngularDistanceComputer angDistComputer) {
        for (final double[] vertex : vpc.vertices(hash, ALL_CARDINALS).values()) {
          final double vLon = vertex[LON_INDEX];
          final double vLat = vertex[LAT_INDEX];
          final double dConeCell = angDistComputer.haversineDistInRad(vLon - coneCenterLon,
              vLat - coneCenterLat, cosCenterLat, cos(vLat));
          if (dConeCell > rRad) {
            return false;
          }
        }
        return true;
      }
    };
    public abstract boolean isOk(double dConeCell, double rRad, VerticesAndPathComputer vpc,
        long hash, double coneCenterLon, double coneCenterLat, double cosCenterLat,
        AngularDistanceComputer angDistComputer);
  }


  public static final class GrowableLongArray {
    public static final Logger LOGGER = Logger.getLogger( NestedSmallCellApproxedMethod.class.getPackage().getName() );
    private long[] array;
    private int cursor;
    public GrowableLongArray(int capacity) {
      assert capacity > 1; // else capacity + (capacity)/2 always returns 1
      this.array = new long[capacity];
      this.cursor = 0;
    }
    public long[] getArray() { return this.array; }
    public int getCursor() { return this.cursor; }
    public final void add(long value) {
      // Mark Taylor suggested to catch the ArrayIndexOutOfBoundsException (I like the idea because 
      // it resorts on the built-in bound check, so we do not have to add an extra test), see 
      //   https://github.com/cds-astro/cds-healpix-java/issues/15
      // I wonder why an explicit test is used in the Java API ArrayList, see e.g.
      //   https://hg.openjdk.java.net/jdk8/jdk8/jdk/file/tip/src/share/classes/java/util/ArrayList.java
      // I guess that it is for better performances when frequent need to make the array grow
      // (but in our case the operation is supposed to be infrequent).
      // Is the compiler/jit smart enough not to make the test (test + bound check) twice?
      // I so far leave as it is: ideally I should have checked performances with both solutions
      // (difference probably negligible, but to be checked!).
      if (this.cursor == this.array.length) {
        // On purpose spurious message (may be better to use an external logger).
        LOGGER.warning("Had to grow moc size! Investigate to find a better estimate!");
        // New size = 1.5 * old size (same code as in Java API ArrayList)
        this.array = Arrays.copyOf(this.array, this.array.length + (this.array.length >> 1));
      }
      this.array[this.cursor++] = value;
    }
  }

  // Because we do not want HealpixNestedLonLatComputer to publicly implement SettableHashParts
  private final SettableHashParts hashPartsProxy = new SettableHashParts() {
    @Override public int baseCellHash() { return baseCellHash; }
    @Override public int iInBaseCell()  { return iInBaseCell; }
    @Override public int jInBaseCell()  { return jInBaseCell; }
    @Override public void setBaseCellHash(int baseCelHash) { baseCellHash = baseCelHash; }
    @Override public void setIInBaseCell(int iInBaseCel)   { iInBaseCell = iInBaseCel; }
    @Override public void setJInBaseCell(int jInBaseCel)   { jInBaseCell = jInBaseCel; }
  };

  public NestedSmallCellApproxedMethod(final int startingDepth, final int deeperDepth, final double radiusRad) {
    assert radiusRad > 0;
    this.startingDepth = startingDepth;
    this.deeperDepth = deeperDepth;
    this.deltaDepthMax = this.deeperDepth -  this.startingDepth;
    this.rRad = radiusRad;
    this.hcc = new VerticesAndPathComputer[this.deltaDepthMax + 1];
    if (this.startingDepth == -1) {
      this.hcc[0] = null;
      this.hComputerStartingDepth = null;
      this.neigSelector = null;
    } else {
      final HealpixNested hn = Healpix.getNested(this.startingDepth);
      this.hcc[0] = hn.newVerticesAndPathComputer();
      this.hComputerStartingDepth = hn.newHashComputer();
      this.neigSelector = hn.newNeighbourSelector();
    }
    this.hnDeeperDepth = Healpix.getNested(this.startingDepth + this.deltaDepthMax);
    this.hComputerDeeperDepth =  this.hnDeeperDepth.newHashComputer();
    for (int i = 1; i <= this.deltaDepthMax; i++) {
      this.hcc[i] = Healpix.getNested(this.startingDepth + i).newVerticesAndPathComputer();
    }
    this.angDistComputer = AngularDistanceComputer.getComputer(this.rRad);
    this.neigList = new FlatHashList(-1, 9); // We don't care about the depth, internal usage only
  }

  @Override
  public double getRadius() {
    return this.rRad;
  }

  @Override
  public HealpixNestedFixedRadiusConeComputer newComputer() {
    return new NestedSmallCellApproxedMethod(this.startingDepth, this.deeperDepth, this.rRad);
  }

  @Override
  public HealpixNestedBMOC overlappingCells(double coneCenterLonRad, final double coneCenterLatRad) {
    return overlapping(coneCenterLonRad, coneCenterLatRad, Mode.OVERLAPPING_CELLS);
  }

  @Override
  public HealpixNestedBMOC overlappingCenters(double coneCenterLonRad, double coneCenterLatRad) {
    return overlapping(coneCenterLonRad, coneCenterLatRad, Mode.OVERLAPPING_CENTERS);
  }

  @Override
  public HealpixNestedBMOC overlappingCells(double coneCenterLonRad, double coneCenterLatRad,
      ReturnedCells returnedCells) {
    switch(returnedCells) {
    case FULLY_IN:
      return overlapping(coneCenterLonRad, coneCenterLatRad, Mode.FULLY_IN);
    case OVERLAPPING:
      return overlappingCells(coneCenterLonRad, coneCenterLatRad);
    case CENTER_IN:
      return overlappingCenters(coneCenterLonRad, coneCenterLatRad);
    default:
      throw new Error("Type " + returnedCells + " not implemented!");
    }
  }

  public HealpixNestedBMOC overlapping(double coneCenterLonRad, double coneCenterLatRad, final Mode mode) {
    // Pre-compute constants
    final double cosConeCenterLat = cos(coneCenterLatRad);
    coneCenterLonRad = normalizeLon(coneCenterLonRad);
    assert 0 <= coneCenterLonRad && coneCenterLonRad <= TWO_PI;
    final GrowableLongArray mocElems = new GrowableLongArray(nMocCellInConeUpperBound());
    int mocSize = 0;
    if (this.startingDepth == -1) {
      for (int h = 0; h < 12; h++) {
        buildMocRecursively(mocElems, 1, h, coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, mode);
      }
    } else {
      // Compute hash of the cell containing the cone center
      final long centerHash = this.hComputerStartingDepth.hash(coneCenterLonRad, coneCenterLatRad); // lat is checked here!
      assert -HALF_PI <= coneCenterLatRad && coneCenterLatRad <= HALF_PI;

      this.neigSelector.neighbours(centerHash, this.neigList);
      this.neigList.put(centerHash);
      this.neigList.sortByHashAsc();

      for (int i = 0; i < this.neigList.size(); i++) {
        buildMocRecursively(mocElems, 0, this.neigList.get(i), coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, mode);
      }
    }
    return HealpixNestedBMOC.createPacking(this.deeperDepth, mocElems.array, mocElems.cursor);
  }

  private final GrowableLongArray buildMocRecursively(final GrowableLongArray moc, int deltaDepth, long hash,
      final double coneCenterLon, final double coneCenterLat, final double cosCenterLat, final Mode mode) {
    final int depth = this.startingDepth + deltaDepth;
    assert this.hcc[deltaDepth].depth() == depth;
    final VerticesAndPathComputer vpc = this.hcc[deltaDepth];
    final double[] center = vpc.center(hash);
    final double cellCenterLon = center[LON_INDEX];
    final double cellCenterLat = center[LAT_INDEX];
    final double dConeCell = this.angDistComputer.haversineDistInRad(cellCenterLon - coneCenterLon,
        cellCenterLat - coneCenterLat, cosCenterLat, cos(cellCenterLat));
    final double rCircumCircle = Healpix.getLargestCenterToCellVertexDistance(
        cellCenterLon, cellCenterLat, depth);
    if (isCellFullyInCone(rRad, rCircumCircle, dConeCell)) {
      moc.add(buildValue(this.startingDepth + deltaDepth, hash, true, this.deeperDepth));
    } else if (isCellOverlapingCone(rRad, rCircumCircle, dConeCell)) {
      if (deltaDepth == this.deltaDepthMax) {
        if (mode.isOk(dConeCell, this.rRad, vpc, hash, coneCenterLon, coneCenterLat, cosCenterLat, this.angDistComputer)) { 
          moc.add(buildValue(this.deeperDepth, hash, false, this.deeperDepth));
        }
      } else {
        hash <<= 2;
        deltaDepth++;
        buildMocRecursively(moc, deltaDepth,   hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
        buildMocRecursively(moc, deltaDepth, ++hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
        buildMocRecursively(moc, deltaDepth, ++hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
        buildMocRecursively(moc, deltaDepth, ++hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
      }
    } // else d > r + c => cell fully out of the cone
    return moc;
  }

  private final int nMocCellInConeUpperBound() {
    // OLD UPPER BOUND: fails in rare circumstances, e.g. depth = 14, radius = 0.001, alpha = 0.002 ,delta = -1.3;
    // cell_area = 4 * pi / ncell = 4 * pi / (3 * 4 * nside^2) = pi / (3 * nside^2) =  pi * r^2
    // cell_radius = r = 1 / (sqrt(3) * nside)
    // As a very simple and naive rule, we take 6x the number of cells needed to cover
    // the cone external annulus
    // Annulus area = 4 pi ((R + r)^2 - R^2) = 4 pi (r^2 + 2rR)
    // N cells = 4 pi (r^2 + 2rR) / 4 pi r^2 = 1 + 2 R/r = 1 + 2 * sqrt(3) * nside * R
    //final double twiceSqrt3 = 2 * 1.73205080756887729352;
    // return 6 * (1 +  (int) (this.hnDeeperDepth.nside * twiceSqrt3 * this.rRad + 0.99));

    // NEW UPPER BOUND: supposedly more robust (and faster to compute)
    // At lower resolution, max 9 cells (partially) overlapped by a cone
    // (except at depth 0, but if more that 9 cells are overlapped, I expect cell(s) to be fully included in the cone) 
    // => grid nside max = 3 * (2^DeltaDepth)
    // For each row (or col) of the higher resolution grid, at the border of the cone, in the worst case: 
    // - 6 = x2 (both sides) x (1 cell overlapping externally + 1 cell overlapping internally + 1 not-merged internal cell)
    // Then at resolution max - 1, #row(or cols) / 2, in the worst case
    // - few chances (?) to have more than 2 consecutive unmerged cells on each side: 4 = 2 * 2 
    // By recursivity, we get:
    // - 6 * 2^DeltaDepth + 4 * 2^(DeltaDepth - 1) + 4 * 2^(DeltaDepth - 2) + ...
    // - 6 * (2^DeltaDepth + 2 * 2^(DeltaDepth + 2^(DeltaDepth + ... )
    // 2^DeltaDepth * (6 + 2 + 1 + 1/2 + 1/4 + ...)
    // We take 12 instead of 10
    return 12 << this.deltaDepthMax;
  }

  private static final boolean isCellFullyInCone(final double coneRadius,
      final double cellCircumCircleRadius, final double coneCenterToCellCenterDistance) {
    return coneCenterToCellCenterDistance <= coneRadius - cellCircumCircleRadius; 
  }

  private static final boolean isCellOverlapingCone(final double coneRadius,
      final double cellCircumCircleRadius, final double coneCenterToCellCenterDistance) {
    return coneCenterToCellCenterDistance < coneRadius + cellCircumCircleRadius; 
  }


  private final int ringIndex(final int deltaDepth, final long hash) {
    this.hnDeeperDepth.decodeRegularHash(hash, this.hashPartsProxy);
    return ringIndex(this.hnDeeperDepth, this.baseCellHash, this.iInBaseCell, this.jInBaseCell);
  }


  private static final int ringIndex(final HealpixNested hn,
      long baseCellHash, int iInBasePixel, int jInBasePixel) {
    final int h = iInBasePixel + jInBasePixel;
    final long jBasePixel = hn.dividedBy4Quotient(baseCellHash);
    return (int) (hn.nsideTime(jBasePixel + 2) - (h + 2));
  }

}
