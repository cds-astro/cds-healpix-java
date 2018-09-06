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
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
// import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_I;
import static cds.healpix.common.math.HackersDelight.toBits;

import java.util.EnumMap;
import java.util.EnumSet;
import java.util.Map;

import cds.healpix.CompassPoint.Cardinal;

final class HealpixNestedVerticesAndPathComputer implements VerticesAndPathComputer {

  // private static final byte[] BASE_CELL_OFFSET_X = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};
  // private static final byte[] BASE_CELL_OFFSET_Y = {1, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1};

  private final HealpixNested h;
  private final HealpixUnprojector unprojector = new HealpixUnprojector();
  private final int nsideMinus1;

  private int baseCellHash;             // Hash value of the base cell (i.e. the depth 0 cell)
  private int iInBaseCell, jInBaseCell; // South-East / South-West coordinates inside a base cell

  private int baseCellOffsetX, baseCellOffsetY; //

  private int xInt, yInt;              // Cell center coordinates in the projection plane
  private double x, y;

  // Because we do not want HealpixNestedLonLatComputer to publicly implement SettableHashParts
  private final SettableHashParts hashPartsProxy = new SettableHashParts() {
    @Override public int baseCellHash() { return baseCellHash; }
    @Override public int iInBaseCell()  { return iInBaseCell; }
    @Override public int jInBaseCell()  { return jInBaseCell; }
    @Override public void setBaseCellHash(int baseCelHash) { baseCellHash = baseCelHash; }
    @Override public void setIInBaseCell(int iInBaseCel)   { iInBaseCell = iInBaseCel; }
    @Override public void setJInBaseCell(int jInBaseCel)   { jInBaseCell = jInBaseCel; }
  };

  HealpixNestedVerticesAndPathComputer(final HealpixNested healpixNested) {
    this.h = healpixNested;
    this.nsideMinus1 = this.h.nsideRemainderMask; // = nside - 1
  }

  @Override
  public int depth() {
    return h.depth;
  }

  @Override
  public double[] center(final long hash) {
    final double[] resultLonLat = new double[2];
    center(hash, resultLonLat);
    return resultLonLat;
  }

  @Override
  public void center(final long hash, final double[] resultLonLat) {
    centerOfPojectedCell(hash, resultLonLat);
    unprojector.unproject(resultLonLat[X_INDEX], resultLonLat[Y_INDEX], resultLonLat);
  }

  @Override
  public double[] vertex(final long hash, final Cardinal cardinalPoint) {
    final double[] resultLonLat = new double[2];
    vertex(hash, cardinalPoint, resultLonLat);
    return resultLonLat;
  }

  @Override
  public void vertex(final long hash, final Cardinal cardinalPoint, final double[] resultLonLat) {
    centerOfPojectedCell(hash, resultLonLat);
    vertexLonLat(resultLonLat, cardinalPoint, resultLonLat);
  }

  @Override
  public EnumMap<Cardinal, double[]> vertices(final long hash, final EnumSet<Cardinal> cardinalPoints) {
    final EnumMap<Cardinal, double[]> verticesMap = new EnumMap<Cardinal, double[]>(Cardinal.class);
    final double[] centerXY = centerOfPojectedCell(hash);
    for (final Cardinal c : cardinalPoints) {
      final double[] resultLonLat = new double[2];
      vertexLonLat(centerXY, c, resultLonLat);
      verticesMap.put(c, resultLonLat);
    }
    return verticesMap;
  }

  @Override
  public void vertices(final long hash, final EnumMap<Cardinal, double[]> cardinalPoints) {
    final double[] centerXY = centerOfPojectedCell(hash);
    for (final Map.Entry<Cardinal, double[]> e : cardinalPoints.entrySet()) {
      vertexLonLat(centerXY, e.getKey(), e.getValue());
    }
  }

  @Override
  public double[][] pathAlongCellSide(final long hash,
      final Cardinal fromVertex, final Cardinal toVertex, 
      boolean isToVertexIncluded, int nSegments) {
    final int resultSize = isToVertexIncluded ? nSegments + 1 : nSegments;
    final double[][] pathPoints = new double[resultSize][];
    pathAlongCellSide(hash, fromVertex, toVertex, isToVertexIncluded, nSegments, pathPoints);
    return pathPoints;
  }

  @Override
  public void pathAlongCellSide(final long hash,
      final Cardinal fromVertex, final Cardinal toVertex, 
      boolean isToVertexIncluded, int nSegments, final double[][] pathPoints) {
    final double[] centerXY = centerOfPojectedCell(hash);
    pathAlongCellSide(centerXY, fromVertex, toVertex, isToVertexIncluded, nSegments, pathPoints, 0);
  }

  @Override
  public double[][] pathAlongCellEdge(final long hash, Cardinal startingVertex,
      final boolean clockwiseDirection, final int nSegmentsBySide) {
    final double[][] pathPoints = new double[nSegmentsBySide << 2][];
    pathAlongCellEdge(hash, startingVertex, clockwiseDirection, nSegmentsBySide, pathPoints);
    return pathPoints;
  }

  @Override
  public void pathAlongCellEdge(final long hash, final Cardinal startingVertex,
      final boolean clockwiseDirection,final int nSegmentsBySide, final double[][] pathPoints) {
    // Compute center
    final double[] centerXY = centerOfPojectedCell(hash);
    // Unrolled loop over successive sides
    // - compute vertex sequence
    Cardinal vertex1 = startingVertex, vertex2, vertex3, vertex4;
    if (clockwiseDirection) {
      vertex2 = startingVertex.nextClockwise();
      vertex3  = vertex2.nextClockwise();
      vertex4 = vertex3.nextClockwise();
    } else {
      vertex2 = startingVertex.nextCounterClockwise();
      vertex3  = vertex2.nextCounterClockwise();
      vertex4 = vertex3.nextCounterClockwise();
    }
    // - make the five sides
    pathAlongCellSide(centerXY, vertex1, vertex2, false, nSegmentsBySide, pathPoints, 0);
    pathAlongCellSide(centerXY, vertex2, vertex3, false, nSegmentsBySide, pathPoints, nSegmentsBySide);
    pathAlongCellSide(centerXY, vertex3, vertex4, false, nSegmentsBySide, pathPoints, nSegmentsBySide << 1);
    pathAlongCellSide(centerXY, vertex4, vertex1, false, nSegmentsBySide, pathPoints, (nSegmentsBySide << 2) - nSegmentsBySide);
  }
  
  private void pathAlongCellSide(final double[] centerXY,
      final Cardinal fromVertex, final Cardinal toVertex, 
      boolean isToVertexIncluded, int nSegments, final double[][] pathPoints, 
      final int fromPathPointsIndex) {
    final int resultSize = isToVertexIncluded ? nSegments + 1 : nSegments;
    // Compute starting point offsets
    final double fromOffsetX = fromVertex.timeXOffset(h.oneOverNside);
    final double fromOffsetY = fromVertex.timeYOffset(h.oneOverNside);
    // Compute stepX and stepY
    final double stepX = (toVertex.timeXOffset(h.oneOverNside) - fromOffsetX) / nSegments;
    final double stepY = (toVertex.timeYOffset(h.oneOverNside) - fromOffsetY) / nSegments;
    // Compute intermediary vertices
    for (int i = 0; i < resultSize; i++) {
      final double[] pathPoint = new double[2];
      final double x = centerXY[X_INDEX] + fromOffsetX + i * stepX;
      final double y = centerXY[Y_INDEX] + fromOffsetY + i * stepY;
      unprojector.unproject(x < 0 ? x + 8: x, y, pathPoint);
      pathPoints[i + fromPathPointsIndex] = pathPoint;
    }
  }

  private void vertexLonLat(final double[] centerXY, final Cardinal vertexDirection,
      final double[] resultVertexLonLat) {
    final double x = centerXY[X_INDEX] + vertexDirection.timeXOffset(h.oneOverNside);
    final double y = centerXY[Y_INDEX] + vertexDirection.timeYOffset(h.oneOverNside);
    unprojector.unproject(x < 0 ? x + 8 : x, y, resultVertexLonLat);
  }


  private double[] centerOfPojectedCell(final long hash) {
    final double[] resultXY = new double[2];
    centerOfPojectedCell(hash, resultXY);
    return resultXY;
  }

  /**
   * 
   * @param hash hash value of the cell we are looking for the coordinates of the center. 
   * @param result the coordinates of the center of the cell in the projection plane,
   *               X in [-pi, pi] and Y in [-pi/2, pi/2].
   */
  private void centerOfPojectedCell(final long hash, final double[] resultXY) {
    checkRegularHash(hash);
    decodeRegularHash(hash);
    rotate45andScale2();
    shiftFromSmallCellCenterToBaseCellCenter();
    scaleToProjectionDividingByNside();
    computeBaseCellCenterOffsetsIn8x3Grid();
    applyBaseCellCenterOffsets();
    setResultIn(resultXY);
  }

  private void checkRegularHash(long hash) {
    if (hash < 0 || h.nHash < hash) {
      throw new IllegalArgumentException("Wrong hash value. Expected value in [0, " + h.nHash +"[;"
          + " Actual value: " + hash);
    }
  }

  private void decodeRegularHash(final long hash) {
    h.decodeRegularHash(hash, hashPartsProxy);
    assert 0 <= baseCellHash && baseCellHash < 12;
    assert 0 <= iInBaseCell && iInBaseCell < h.nside;
    assert 0 <= jInBaseCell && jInBaseCell < h.nside;
  }

  private void rotate45andScale2() {
    // Rotation +45 deg and scale of factor2lBaseCell
    xInt = iInBaseCell - jInBaseCell;    assert -h.nside <  xInt && xInt < h.nside;
    yInt = iInBaseCell + jInBaseCell;    assert        0 <= yInt && yInt <= 2 * (h.nside - 1);
  }

  private void shiftFromSmallCellCenterToBaseCellCenter() {
    yInt -= nsideMinus1;    assert -h.nside < yInt && yInt < h.nside;
  }

  private void scaleToProjectionDividingByNside() {
    x = xInt * h.oneOverNside;
    y = yInt * h.oneOverNside;
  }

  private void computeBaseCellCenterOffsetsIn8x3Grid() {
    // baseCellOffsetX = BASE_CELL_OFFSET_X[baseCellHash];
    // baseCellOffsetY = BASE_CELL_OFFSET_Y[baseCellHash];
    baseCellOffsetY = divBy4Quotient(baseCellHash);  assert 0 <= baseCellOffsetY && baseCellOffsetY <= 2;
    baseCellOffsetY = 1 - baseCellOffsetY;           assert -1 <= baseCellOffsetY && baseCellOffsetY <= 1;

    baseCellOffsetX = divBy4Reminder(baseCellHash);  assert 0 <= baseCellOffsetX && baseCellOffsetX <= 3;   
    baseCellOffsetX <<= 1;                           assert baseCellOffsetX == 0 || baseCellOffsetX == 2
                                                         || baseCellOffsetX == 4 || baseCellOffsetX == 6;
   // +1 if the base cell is not equatorial
    baseCellOffsetX |= baseCellOffsetY & 1;          assert 0 <= baseCellOffsetX && baseCellOffsetX <= 7;
   // +8 if the base cell is number 4 and xInt is negative
   // baseCellOffsetX |= ((xInt & SIGN_BIT_MASK_I) >>> (24 | baseCellHash)) & 8;
  }


  private static int divBy4Quotient(int val) {
    return val >> 2;
  }
  private static int divBy4Reminder(int val) {
    return val & 3;
  }

  private void applyBaseCellCenterOffsets() {
    x += baseCellOffsetX;
    y += baseCellOffsetY;
    // If x < 0, then x += 8; (happens only in case of base cell 4)
    x += (toBits(x) & SIGN_BIT_MASK_L) >>> 60;
  }

  private void setResultIn(final double[] resultXY) {
    resultXY[X_INDEX] = x;
    resultXY[Y_INDEX] = y;
  }
}
