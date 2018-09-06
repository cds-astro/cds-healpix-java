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
import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;
import static cds.healpix.common.math.Math.cos;

final class NestedLargeCellApproxedMethod implements HealpixNestedFixedRadiusConeComputer {

  private final int startingDepth;
  private final int hashDepth;
  private final HealpixNested hn;
  private final HashComputer hc;
  private final NeighbourSelector neiSelect;
  private final VerticesAndPathComputer vpc;
  private final int twiceDeltaOrder;
  private final double coneRadiusRad;
  
  private final AngularDistanceComputer angDistComputer;
  
  private final long[] result = new long[9];
  private final double[] cellCenter = new double[2];
  private final NeighbourList neigList;
  
  public NestedLargeCellApproxedMethod(final int startingDepth, final int outputDepth, final double radiusRad) {
    this.startingDepth = startingDepth;
    this.hashDepth = outputDepth;
    if (startingDepth < this.hashDepth) {
      throw new IllegalArgumentException("Choosen radius too large for the wanted depth!");
    }
    this.hn = Healpix.getNested(startingDepth);
    this.hc = this.hn.newHashComputer();
    this.vpc = this.hn.newVerticesAndPathComputer();
    this.neiSelect = this.hn.newNeighbourSelector();
    this.twiceDeltaOrder = (startingDepth - hashDepth) << 1;
    this.coneRadiusRad = radiusRad;
    this.angDistComputer = AngularDistanceComputer.getComputer(this.coneRadiusRad);
    this.neigList = new NeighbourList(-1); // We don't care, internal usage only.
  }
  
  @Override
  public double getRadius() {
    return this.coneRadiusRad;
  }

  @Override
  public HealpixNestedFixedRadiusConeComputer newComputer() {
    return new NestedLargeCellApproxedMethod(this.startingDepth, this.hashDepth, this.coneRadiusRad);
  }
  
  @Override
  public HealpixNestedBMOC overlappingCells(double coneCenterLon, final double coneCenterLat) {
    final double cosConeCenterLat = cos(coneCenterLat);
    
    final long centerHash = this.hc.hash(coneCenterLon, coneCenterLat);
    this.neiSelect.neighbours(centerHash, this.neigList);
    final long ch = centerHash >>> this.twiceDeltaOrder;
    final long ch4moc = buildValue(this.hashDepth, ch, false, this.hashDepth);
    result[0] = ch4moc;
    int ir = 1;
    for (int i = 0; i < this.neigList.size(); i++) {
      final long h = this.neigList.get(i);
      final long hmm = h >>> this.twiceDeltaOrder;
      final long h4moc = buildValue(this.hashDepth, hmm, false, this.hashDepth);
      if (h4moc != ch4moc && isNotIn(h4moc, result, ir)) {
        this.vpc.center(h, cellCenter);
        final double cellCenterLon = cellCenter[LON_INDEX];
        final double cellCenterLat = cellCenter[LAT_INDEX];
        final double dConeCell = this.angDistComputer.haversineDistInRad(cellCenterLon - coneCenterLon,
            cellCenterLat - coneCenterLat, cosConeCenterLat, cos(cellCenterLat));
        if (dConeCell <= this.coneRadiusRad) {
          result[ir++] = h4moc; 
        }
      }
    }
    return HealpixNestedBMOC.createPacking(this.hashDepth, result, ir);
  }
  
  @Override
  public HealpixNestedBMOC overlappingCenters(double coneCenterLon, double coneCenterLat) {
final double cosConeCenterLat = cos(coneCenterLat);
    
    final long centerHash = this.hc.hash(coneCenterLon, coneCenterLat);
    this.neiSelect.neighbours(centerHash, this.neigList);
    final long ch = centerHash >>> this.twiceDeltaOrder;
    final long ch4moc = buildValue(this.hashDepth, ch, false, this.hashDepth);
    result[0] = ch4moc;
    int ir = 1;
    for (int i = 0; i < this.neigList.size(); i++) {
      final long h = this.neigList.get(i);
      final long hmm = h >>> this.twiceDeltaOrder;
      final long h4moc = buildValue(this.hashDepth, hmm, false, this.hashDepth);
      if (h4moc != ch4moc && isNotIn(h4moc, result, ir)) {
        this.vpc.center(h, cellCenter);
        final double cellCenterLon = cellCenter[LON_INDEX];
        final double cellCenterLat = cellCenter[LAT_INDEX];
        final double rCircumCircle = Healpix.getLargestCenterToCellVertexDistance(
          cellCenterLon, cellCenterLat, this.startingDepth);
        final double dConeCell = this.angDistComputer.haversineDistInRad(cellCenterLon - coneCenterLon,
            cellCenterLat - coneCenterLat, cosConeCenterLat, cos(cellCenterLat));
        if (isCellOverlapingCone(this.coneRadiusRad, rCircumCircle, dConeCell)) {
          result[ir++] = h4moc; 
        }
      }
    }
    return HealpixNestedBMOC.createPacking(this.hashDepth, result, ir);
  }
  
  
  private static final boolean isNotIn(final long h, final long[] ah, final long alength) {
    for (int i = 1; i < alength; i++) {
      if (h == ah[i]) {
        return false;
      }
    }
    return true;
  }

  private static final boolean isCellOverlapingCone(final double coneRadius,
      final double cellCircumCircleRadius, final double coneCenterToCellCenterDistance) {
   return coneCenterToCellCenterDistance < coneRadius + cellCircumCircleRadius; 
  }
  
}
