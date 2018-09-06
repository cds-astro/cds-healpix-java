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

import static cds.healpix.HealpixNestedBMOC.buildValue;
import static cds.healpix.NestedSmallCell.normalizeLon;
import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.cos;

final class NestedSmallCellApproxedMethod implements HealpixNestedFixedRadiusConeComputer {

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
      public boolean isOk(final double dConeCell, final double rRad) {
        return true;
      }
    },
    OVERLAPPING_CENTERS() {
      @Override
      public boolean isOk(final double dConeCell, final double rRad) {
        return dConeCell <= rRad;
      }
    };
    public abstract boolean isOk(double dConeCell, final double rRad);
  }
  
  
  //Because we do not want HealpixNestedLonLatComputer to publicly implement SettableHashParts
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
  
  public HealpixNestedBMOC overlapping(double coneCenterLonRad, double coneCenterLatRad, final Mode mode) {
    // Pre-compute constants
    final double cosConeCenterLat = cos(coneCenterLatRad);
    coneCenterLonRad = normalizeLon(coneCenterLonRad);
    assert 0 <= coneCenterLonRad && coneCenterLonRad <= TWO_PI;
    // Store required space in the MOC
    final long[] mocElems = new long[mocSizeUpperLimit(coneCenterLonRad, coneCenterLatRad)];
    int mocSize = 0;
    if (this.startingDepth == -1) {
      for (int h = 0; h < 12; h++) {
        mocSize = buildMocRecursively(mocElems, mocSize, 1, h, coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, mode);
      }
    } else {
      // Compute hash of the cell containing the cone center
      final long centerHash = this.hComputerStartingDepth.hash(coneCenterLonRad, coneCenterLatRad); // lat is checked here!
      assert -HALF_PI <= coneCenterLatRad && coneCenterLatRad <= HALF_PI;
     
      this.neigSelector.neighbours(centerHash, this.neigList);
      this.neigList.put(centerHash);
      this.neigList.sortByHashAsc();

      for (int i = 0; i < this.neigList.size(); i++) {
        mocSize = buildMocRecursively(mocElems, mocSize, 0, this.neigList.get(i), coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, mode);
      }
    }
    return HealpixNestedBMOC.createPacking(this.deeperDepth, mocElems, mocSize);
  }
  
  private final int buildMocRecursively(final long[] moc, int mocLength, int deltaDepth, long hash,
      final double coneCenterLon, final double coneCenterLat, final double cosCenterLat, final Mode mode) {
    final int depth = this.startingDepth + deltaDepth;
    assert this.hcc[deltaDepth].depth() == depth;
    final double[] center = this.hcc[deltaDepth].center(hash);
    final double cellCenterLon = center[LON_INDEX];
    final double cellCenterLat = center[LAT_INDEX];
    final double dConeCell = this.angDistComputer.haversineDistInRad(cellCenterLon - coneCenterLon,
        cellCenterLat - coneCenterLat, cosCenterLat, cos(cellCenterLat));
    final double rCircumCircle = Healpix.getLargestCenterToCellVertexDistance(
        cellCenterLon, cellCenterLat, depth);
    if (isCellFullyInCone(rRad, rCircumCircle, dConeCell)) {
      moc[mocLength++] = buildValue(this.startingDepth + deltaDepth, hash, true, this.deeperDepth);
    } else if (isCellOverlapingCone(rRad, rCircumCircle, dConeCell) && mode.isOk(dConeCell, this.rRad)) {
      if (deltaDepth == this.deltaDepthMax) {
        moc[mocLength++] = buildValue(this.deeperDepth, hash, false, this.deeperDepth);
      } else {
        hash <<= 2;
        deltaDepth++;
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth,   hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, coneCenterLon, coneCenterLat, cosCenterLat, mode);
      }
    } // else d > r + c => cell fully out of the cone
    return mocLength;
  }
  
  private final int mocSizeUpperLimit(final double coneCenterLonRad, final double coneCenterLatRad) {
    // Count the number of rings in the other sie of a pole (if the cone contain a pole),
    // the smallest possible ring index...
    int nRingsOtherSideOfPole = 0;
    int smallestCornerRingIndex = -1;  
    double latMax = coneCenterLatRad + this.rRad;
    if (latMax > HALF_PI) {
      nRingsOtherSideOfPole = ringIndex(this.deltaDepthMax, this.hComputerDeeperDepth.hash(coneCenterLonRad + PI, HALF_PI - latMax)) + 1;
      latMax = HALF_PI;
    } else {
      smallestCornerRingIndex // min = -1 = North polevpComputers
      = Math.max(smallestCornerRingIndex,
          ringIndex(this.deltaDepthMax, this.hComputerDeeperDepth.hash(coneCenterLonRad, latMax)) - 1);
    }
    // ... and the largest possible ring index
    int largestCornerRingIndex = Healpix.nIsolatitudeRings(this.deeperDepth);
    double latMin = coneCenterLatRad - this.rRad;
    if (latMin < -HALF_PI) {
      nRingsOtherSideOfPole = largestCornerRingIndex - ringIndex(this.deltaDepthMax, this.hComputerDeeperDepth.hash(coneCenterLonRad + PI, -HALF_PI - latMin));
      latMin = -HALF_PI;
    } else {
      largestCornerRingIndex = Math.min(largestCornerRingIndex, ringIndex(this.deltaDepthMax, this.hComputerDeeperDepth.hash(coneCenterLonRad, latMin)) + 1); 
    }
    // Deduce the number of rings...
    final int nRings = largestCornerRingIndex - smallestCornerRingIndex + 1;
    // ... and the upper limit on the MOC size
    return nElemMax(nRings, nRingsOtherSideOfPole);
  }
  
  private static final boolean isCellFullyInCone(final double coneRadius,
      final double cellCircumCircleRadius, final double coneCenterToCellCenterDistance) {
   return coneCenterToCellCenterDistance <= coneRadius - cellCircumCircleRadius; 
  }
  
  private static final boolean isCellOverlapingCone(final double coneRadius,
      final double cellCircumCircleRadius, final double coneCenterToCellCenterDistance) {
   return coneCenterToCellCenterDistance < coneRadius + cellCircumCircleRadius; 
  }
  
  private static final int nElemMax(final int nRings, final int nRingOtherSideOfThePole) {
    //   nRing * 2      (smallest cells at the extremity of each ring)
    // + nRing / 2 * 2  (cells of depth deppestDepth - 1 a thte extremity of each ring)
    // + nRing / 4 * 2
    // + ...
    // + nRing / DeltaDepth * 2
    // = nRing * (2 + 1 + 1/2 + 1/4 + ...)
    // and the serie 1/2 + 1/4 + ... = 1
    // ~= nRing * (2 + 1 + 1) = 4  * nRing
    // Constant should be 4, we use 8 to be conservative
    return 8 * (nRings + nRingOtherSideOfThePole);
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
