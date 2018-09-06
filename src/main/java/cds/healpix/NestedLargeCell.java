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

import java.util.EnumMap;
import java.util.EnumSet;

import cds.healpix.common.math.HackersDelight;

import cds.healpix.CompassPoint.MainWind;
import cds.healpix.CompassPoint.Cardinal;

import static cds.healpix.CompassPoint.Cardinal.E;
import static cds.healpix.CompassPoint.Cardinal.N;
import static cds.healpix.CompassPoint.Cardinal.S;
import static cds.healpix.CompassPoint.Cardinal.W;

import static cds.healpix.Healpix.isLatInNorthPolarCap;
import static cds.healpix.Healpix.isLatInSouthPolarCap;

import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;

import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.ONE_OVER_2PI;
import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.abs;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.asin;
import static cds.healpix.common.math.Math.isFinite;

/**
 * Special case in which the radius of the cone is smaller than the smallest possible distance
 * to the edge of a neighbouring cell at the requested depth => the cone overlap with maximum 9
 * cells (the cell containing the center of the cone plus the (max) 8 neighbours.
 * If the cone contains one corner
 *   => we returns the five cells sharing this same corner
 * If rCone + distance(coneCenter, pixelCenter) < smalllestDistanceFromCenterToEdge
 *   => only the cell containing the cone center is to be returned
 * Else, we look at the quadran (NE, NW, SW or SE) of the center of the cone insquareOfSinOfHalfRadius the cell and
 *   look at the possible intersection with the border of the quadran
 *   
 * In the case of small angles, we can compute the intersection between an ellipse and a strait line
 * NE: y = -x + b
 * NW: y =  x + b
 * SW: y = -x + b
 * SE: y =  x + b
 * 
 * b can be computed from a corner: yc = xc + b => b = (yc - xc)
 * 
 * To do so, we compute the equation of the strait line in the space in which the ellipse is a
 * circle of radius one.
 * 1 - put the center of the frame at the center of the cone => x = (x - xc) and y = (y - yc)
 *   => y = x + b => y = x + (b+yc-xc)
 * 2 - if var/cov matrix not diagonal, diagonalize and deduce the rotation matrix M
 * 
 * 3 - reduce x by sigx and y by sigy
 *   => y = x + b, y' = y/sigy, x' = x/sigx
 *   => y = y'*sigy = x'*sigx + b
 *   => y' = (sigx * x' + b) / sigy
 *   or
 *   => sigy * y' -sigx * x' - b = 0 (for by+ax+c=0)
 *   
 * 4 - diqtence entre un point et une droite
 *   d = |sigy * y'c - sigx * x'c - b| / sqrt(sigy^2+sigx^2)
 *   but y'c = x'c = 0
 *   => d = b/sqrt()
 * 
 * 
 * OK if cone radius > 20 mas !!!!
 * Else, MAKE THE FLAT SKY APPROXIMATION, USE VERTICES COORDINATES TO COMPUTE INTERSECTIONS!!
 * Orthographic projection: x = cos(lat - latc) * sin(lon - lonc)
 *                          y = sin(lat - latc)
 * With cos(lat - latc) ~= 1
 *      sin(a-b) ~= a - b
 * 
 * @author F.-X. Pineau
 *
 */
final class NestedLargeCell implements HealpixNestedFixedRadiusConeComputer {

  // private static final Logger LOG = LoggerFactory.getLogger("fr.unistra.cds.healpix.cone"); // NestedLargeCell.class.getName()

  private final int nIterMax;
  private final double epsZ;

  private final HealpixNested hn;
  private final double rRad;
  private final double twoSineOfHalfConeRadius;
  private final double squareOfsinOfHalfR;


  private final HashComputer hComputer;
  private final NeighbourSelector neigSelector;
  private final VerticesAndPathComputer vpComputer;

  private final AngularDistanceComputer angDistComputer;

  private final EnumMap<Cardinal, double[]> vertices = new EnumMap<Cardinal, double[]>(Cardinal.class);

  private final EnumSet<MainWind> neigSet = EnumSet.noneOf(MainWind.class);
  private final NeighbourList neigList;

  private final long[] neig = new long[8];
  private int iNeig = 0;
  
  public NestedLargeCell(final HealpixNested hn, final double radiusRad) {
    this(hn, radiusRad, radiusRad * 1e-4, 10);
  }

  protected NestedLargeCell(final HealpixNested hn, final double radiusRad,
      final double eps, final int nIterMax) {
    assert radiusRad > 0;
    this.hn = hn;
    this.rRad = radiusRad;
    this.epsZ = eps;
    this.nIterMax = nIterMax;
    this.angDistComputer = AngularDistanceComputer.getComputer(this.rRad);
    final double sinHalfR = angDistComputer.sin(0.5 * radiusRad); 
    this.twoSineOfHalfConeRadius = 2 * sinHalfR;
    this.squareOfsinOfHalfR = sinHalfR * sinHalfR;

    this.hComputer = this.hn.newHashComputer();
    this.neigSelector = this.hn.newNeighbourSelector();
    this.vpComputer = this.hn.newVerticesAndPathComputer();
    this.neigList = new NeighbourList(-1); // We don't care, internal usage only
    initVerticesMap();
  }

  private void initVerticesMap() {
    for (final Cardinal c : Cardinal.values()) {
      vertices.put(c, new double[2]);
    }
  }

  @Override
  public double getRadius() {
    return this.rRad;
  }

  @Override
  public HealpixNestedBMOC overlappingCenters(double coneCenterLonRad, double coneCenterLatRad) {
    final double cosConeCenterLat = cos(coneCenterLatRad);
    // Compute hash of the cell containing the cone center
    final long centerHash = this.hComputer.hash(coneCenterLonRad, coneCenterLatRad); // lat is checked here!
    this.neigSelector.neighbours(centerHash, neigList);
    // Build result based on vertices
    this.neig[0] = centerHash;
    this.iNeig = 1;
    final double[] center = new double[2];
    for (int i = 0; i < this.neigList.size(); i++) {
      final long neigH = this.neigList.get(i);
      this.vpComputer.center(neigH, center);
      final double lat = center[LAT_INDEX];
      final double deltaLon = computeDeltaLon(coneCenterLonRad, center[LON_INDEX]);
      final double deltaLat = coneCenterLatRad - lat;
      if (squareOfSinOfHalfR(deltaLon, deltaLat, cos(lat), cosConeCenterLat) <= squareOfsinOfHalfR) {
        insertSortRmDuplicates(neigH);
      }
    }
    // assert iNeig == this.neigList.size() + 1;
    // Build moc from set set of Hash
    final long[] moc = new long[this.iNeig];
    for (int i = 0; i < this.iNeig; i++) {
      moc[i] = HealpixNestedBMOC.buildValue(this.hn.depth(), neig[i], false, this.hn.depth());
    }
    return HealpixNestedBMOC.createUnsafe(this.hn.depth(), moc);
  }

  
  @Override
  public HealpixNestedBMOC overlappingCells(double coneCenterLonRad, double coneCenterLatRad) {    
    final double cosConeCenterLat = cos(coneCenterLatRad);
    final double sinConeCenterLat = sin(coneCenterLatRad);

// FIRST TEST IF IS IN CIRCLE?? I.E. SMALLEST DISTANCE FROM CENTER TO NEAREST EDGE
    
    coneCenterLonRad = normalizeLon(coneCenterLonRad);
    // LOG.debug("Cone center (long, lat): ({}, {}); radius: {} rad; lon after normalization: {}",
    //     coneCenterLonRad, coneCenterLatRad, this.rRad, coneCenterLonRad);
    assert 0 <= coneCenterLonRad && coneCenterLonRad <= TWO_PI;

    // Compute hash of the cell containing the cone center
    final long centerHash = this.hComputer.hash(coneCenterLonRad, coneCenterLatRad); // lat is checked here!
    // LOG.debug("Hash of the cone center: {}, depth: {}", centerHash, this.hComputer.depth());
    assert -HALF_PI <= coneCenterLatRad && coneCenterLatRad <= HALF_PI;


    // Get the coordinates of the four vertices
    this.vpComputer.vertices(centerHash, this.vertices);

    final double latMax = coneCenterLatRad + this.rRad;
    final double latMin = coneCenterLatRad - this.rRad;
    // LOG.debug("Lat min: {} deg; Lat max: {} deg.", Math.toDegrees(latMin), Math.toDegrees(latMax));
    int nVerticesInCone = 0;
    boolean testNE = true;
    boolean testNW = true;
    boolean testSE = true;
    boolean testSW = true;
    boolean testN = true;
    boolean testS = true;
    // Test if North vertex is in the cone
    final double northVertexLat = this.vertices.get(N)[LAT_INDEX];
    if (northVertexLat <= latMax) {
      final double deltaLon = computeDeltaLon(coneCenterLonRad, this.vertices.get(N)[LON_INDEX]);
      final double deltaLat = northVertexLat - coneCenterLatRad;
      final double cosNorthVertexLat = cos(northVertexLat);
      if (squareOfSinOfHalfR(deltaLon, deltaLat, cosNorthVertexLat, cosConeCenterLat)
          <= squareOfsinOfHalfR) {
        neigSet.add(MainWind.N);
        neigSet.add(MainWind.NE);
        neigSet.add(MainWind.NW);
        nVerticesInCone++;
        testN = false;
        testNE = false;
        testNW = false;
        // LOG.debug("Add N");
      }
    }
    // Test if South vertex is in the cone
    final double southVertexLat = this.vertices.get(S)[LAT_INDEX];
    if (latMin <= southVertexLat) {
      final double deltaLon = computeDeltaLon(coneCenterLonRad, this.vertices.get(S)[LON_INDEX]);
      final double deltaLat = coneCenterLatRad - southVertexLat;
      final double cosSouthVertexLat = cos(southVertexLat);
      if (squareOfSinOfHalfR(deltaLon, deltaLat, cosSouthVertexLat, cosConeCenterLat)
          <= squareOfsinOfHalfR) {
        neigSet.add(MainWind.S);
        neigSet.add(MainWind.SE);
        neigSet.add(MainWind.SW);
        nVerticesInCone++;
        testS = false;
        testSE = false;
        testSW = false;
        // LOG.debug("Add S");
      }
    }
    // Test if East and West vertices are in the cone
    final double eastWestVerticesLat = this.vertices.get(W)[LAT_INDEX];
    assert this.vertices.get(W)[LAT_INDEX] == this.vertices.get(E)[LAT_INDEX];
    if (latMin <= eastWestVerticesLat && eastWestVerticesLat <= latMax) {
      final double cosEastWestVerticesLat = cos(eastWestVerticesLat);
      final double squareOfSinOfHalfDeltaLonMax = squareOfSinOfHalfDeltaLon(
          eastWestVerticesLat - coneCenterLatRad, cosEastWestVerticesLat, cosConeCenterLat);
      // Test East
      if (squareOfSinOfHalfOfDeltaLon(computeDeltaLon(coneCenterLonRad, this.vertices.get(E)[LON_INDEX]))
          <= squareOfSinOfHalfDeltaLonMax) {
        neigSet.add(MainWind.E);
        neigSet.add(MainWind.NE);
        neigSet.add(MainWind.SE);
        nVerticesInCone++;
        testNE = false;
        testSE = false;
        // LOG.debug("Add E");
      }
      // Test West
      if (squareOfSinOfHalfOfDeltaLon(computeDeltaLon(coneCenterLonRad, this.vertices.get(W)[LON_INDEX]))
          <= squareOfSinOfHalfDeltaLonMax) {
        neigSet.add(MainWind.W);
        neigSet.add(MainWind.NW);
        neigSet.add(MainWind.SW);
        nVerticesInCone++;
        testNW = false;
        testSW = false;
     // LOG.debug("Add W");
      }
    }
    // Build result based on vertices
    this.neig[0] = centerHash;
    this.iNeig = 1;
    if (nVerticesInCone > 0) {
      this.neigSelector.neighbours(centerHash, this.neigSet, this.neigList);
      for (int i = 0; i < this.neigList.size(); i++) {
        insertSortRmDuplicates(this.neigList.get(i));
      }
      assert iNeig == this.neigList.size() + 1;
    }
    // Test edge intersections
    if (testNE || testNW || testN) { // NE / NW edges
      if (isLatInNorthPolarCap(northVertexLat)) { // NE/NW edges in the NPC
        edgesNEWinNPC(coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, sinConeCenterLat,
            testNE, testNW, testN);
      } else if (isLatInSouthPolarCap(eastWestVerticesLat)) { // NE/NW edges in the SPC
        edgesNEWinSPC(coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, sinConeCenterLat,
            testNE, testNW);
      } else if (testNE || testNW) { //  NE/NW edges in the EQR
        edgesNEWinEQR(coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, sinConeCenterLat,
            testNE, testNW);
      }
    }
    if (testSE || testSW || testS) { // SE / SW edges
      if (isLatInNorthPolarCap(eastWestVerticesLat)) {   // SE/SW edges in north polar cap
        edgesSEWinNPC(coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, sinConeCenterLat,
            testNE, testNW);
      } else if (isLatInSouthPolarCap(southVertexLat)) { // SE/SW edges in south polar cap
        edgesSEWinSPC(coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, sinConeCenterLat,
            testNE, testNW, testS);
      } else if (testSE || testSW) { // SE/SW edges in the EQR
        edgesSEWinEQR(coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, sinConeCenterLat,
            testNE, testNW);
      }
    }
    // Build moc from set set of Hash
    final long[] moc = new long[this.iNeig];
    for (int i = 0; i < this.iNeig; i++) {
      moc[i] = HealpixNestedBMOC.buildValue(this.hn.depth(), neig[i], false, this.hn.depth());
    }
    return HealpixNestedBMOC.createUnsafe(this.hn.depth(), moc);
    // Get hashes of neighbour cells
    /*this.neigSelector.neighbours(centerHash, this.neigSet, this.neigList);
    // Build moc from the nieghbour list
    long[] moc = new long[this.neigList.size() + 1];
    moc[0] = HealpixNestedMOC.buildValue(this.hn.depth(), centerHash, nVerticesInCone == 4, this.hn.depth());
    int i = 0;
    for (final long hash : this.neigList) {
      moc[++i] = HealpixNestedMOC.buildValue(this.hn.depth(), hash, false, this.hn.depth());
    }
    Arrays.sort(moc);
    return HealpixNestedMOC.createUnsafe(this.hn.depth(), moc);*/
  }

  private final void insertSortRmDuplicates(final long hash) {
    for (int i = 0; i < iNeig; i++) { // No need for a binary search on such a few number of elements
      final long hashi = neig[i];
      if (hash == hashi) { // No duplicates
        return;
      } else if (hash < hashi) { // Insert in the middle of the array
        for (int j = iNeig++; j > i;) {
          neig[j] = neig[--j];
        }
        neig[i] = hash;
        return;
      }
    }
    // Insert a the end of the array
    neig[iNeig++] = hash;
  }

  private final void edgesNEWinEQR(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testNE, final boolean testNW) {
    final double zDxDyEq1 = NewtonMethod.newtonSolveEquatorialZone(
        Math.sin(coneCenterLatRad + rRad * 0.9),
        sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax);
    final double pointLat = asin(zDxDyEq1);
    final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        pointLat - coneCenterLatRad, zDxDyEq1);
    //   LOG.debug("Test N-E/W edges in EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
    //     zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
    final boolean isValidDeltaLon = isFinite(deltaLon);
    if (testNE && isValidDeltaLon) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
    }
    if (testNW && isValidDeltaLon) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
    } 
  }

  private final void edgesSEWinEQR(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testSE, final boolean testSW) {
    final double zDxDyEq1 = NewtonMethod.newtonSolveEquatorialZone(
        Math.sin(coneCenterLatRad - rRad * 0.9),
        sinConeCenterLat, twoSineOfHalfConeRadius, true, epsZ, nIterMax);
    final double pointLat = asin(zDxDyEq1); 
    final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - pointLat, zDxDyEq1);
 // LOG.debug("Test SEW edge in EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
 //     zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
    final boolean isValidDeltaLon = isFinite(deltaLon);
    if (testSE && isValidDeltaLon) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
    }
    if (testSW && isValidDeltaLon) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
    }
  }

  private final void edgesNEWinNPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testNE, final boolean testNW, final boolean testN) {
    // Remark: North pole not included, else (isNotSelectedNE || isNotSelectedNW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    if (testNE || testN) {
      double zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(coneCenterLatRad + rRad * 0.9), true, coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, true, epsZ, nIterMax);
      // TODO: check precisions in case z near from 1 ?
      double pointLat = asin(zDxDyEq1); 
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          pointLat - coneCenterLatRad, zDxDyEq1);
      // LOG.debug("Test NE edge in NPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
      //    zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
      }
      if (coneCenterLonModHalfPi + deltaLon > PI_OVER_FOUR) { // Change the base pixel!
       zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
            Math.sin(coneCenterLatRad + rRad * 0.9), true, coneCenterLonModHalfPi - HALF_PI,
            sinConeCenterLat, twoSineOfHalfConeRadius, true, epsZ, nIterMax);
            /*zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
                Math.sin(coneCenterLatRad - rRad * 0.9), true, coneCenterLonModHalfPi - HALF_PI,
                sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax); // SE!!*/
        pointLat = asin(zDxDyEq1);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            pointLat - coneCenterLatRad, zDxDyEq1);
        // LOG.debug("Test NE edge in NPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
        //    zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
        if (isFinite(deltaLon)) {
          insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
        }
      }
    }
    if (testNW || testN) {
      double zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(coneCenterLatRad + rRad * 0.9), false, coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax);
      double pointLat = asin(zDxDyEq1);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          pointLat - coneCenterLatRad, zDxDyEq1);
      // LOG.debug("Test NW edge in NPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
      //     zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
      }
      if (coneCenterLonModHalfPi - deltaLon < 0) { // Change the base pixel!
        zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
            Math.sin(coneCenterLatRad + rRad * 0.9), false, HALF_PI + coneCenterLonModHalfPi,
            sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax);
        pointLat = asin(zDxDyEq1);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            pointLat - coneCenterLatRad, zDxDyEq1);
        // E.g. in Aladin: draw circle(089.91417, +47.04352, 59')
     // LOG.debug("Test NW edge in NPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
     //     zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
        if (isFinite(deltaLon)) {
          insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
        }
      }
    }
  }

  private final void edgesNEWinSPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testNE, final boolean testNW) {
    // Remark: North pole not included, else (isNotSelectedNE || isNotSelectedNW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    if (testNE) {
      final double zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(-coneCenterLatRad - rRad * 0.9), true, coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax);
      // TODO: check precisions in case z near from 1 ?
      final double pointLat = asin(zDxDyEq1); 
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          pointLat - coneCenterLatRad, zDxDyEq1);
      // LOG.debug("Test NE edge in SPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
      //    zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
      }
    }
    if (testNW) {
      final double zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(-coneCenterLatRad - rRad * 0.9), false, coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, true, epsZ, nIterMax);
      final double pointLat = asin(zDxDyEq1);
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          pointLat - coneCenterLatRad, zDxDyEq1);
   // LOG.debug("Test NW edge in SPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
   //     zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
      }
    }
  }
  
  private final void edgesSEWinNPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testSE, final boolean testSW) {
    // Remark: South pole not included, else (testSE || testSW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    if (testSE) {
      final double zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(coneCenterLatRad - rRad * 0.9), true, coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax);
      // TODO: check precisions in case z near from -1 ?
      final double pointLat = asin(zDxDyEq1); 
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - pointLat, zDxDyEq1);
   // LOG.debug("Test SE edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
   //     zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
      }
    }
    if (testSW) {
      final double zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(coneCenterLatRad - rRad * 0.9), false, coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, true, epsZ, nIterMax);
      final double pointLat = asin(zDxDyEq1); 
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - pointLat, zDxDyEq1);
   // LOG.debug("Test SW edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
   //     zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
      }
    }
  }
  
  private final void edgesSEWinSPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testSE, final boolean testSW, final boolean testS) {
    // Remark: South pole not included, else (testSE || testSW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    if (testSE || testS) {
      double zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(-coneCenterLatRad + rRad * 0.9), true, coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, true, epsZ, nIterMax);
      // TODO: check precisions in case z near from -1 ?
      double pointLat = asin(zDxDyEq1); 
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - pointLat, zDxDyEq1);
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
      }
      if (coneCenterLonModHalfPi + deltaLon > PI_OVER_FOUR) {
        zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
            Math.sin(-coneCenterLatRad + rRad * 0.9), true, coneCenterLonModHalfPi - HALF_PI,
            -sinConeCenterLat, twoSineOfHalfConeRadius, true, epsZ, nIterMax);
        pointLat = asin(zDxDyEq1); 
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            coneCenterLatRad - pointLat, zDxDyEq1);
     // LOG.debug("Test SE edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
     //    zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
        if (isFinite(deltaLon)) {
          insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, pointLat));
        }
      }
    }
    if (testSW || testS) {
      double zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          Math.sin(-coneCenterLatRad + rRad * 0.9), false, coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax);
      double pointLat = asin(zDxDyEq1); 
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - pointLat, zDxDyEq1);
      if (isFinite(deltaLon)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
      }
      if (coneCenterLonModHalfPi - deltaLon < 0) {
        zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
            Math.sin(-coneCenterLatRad + rRad * 0.9), false,  coneCenterLonModHalfPi + HALF_PI,
            -sinConeCenterLat, twoSineOfHalfConeRadius, false, epsZ, nIterMax);
        pointLat = asin(zDxDyEq1); 
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            coneCenterLatRad - pointLat, zDxDyEq1);
     // LOG.debug("Test SW edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg",
     //    zDxDyEq1, Math.toDegrees(pointLat), Math.toDegrees(deltaLon));
        if (isFinite(deltaLon)) {
          insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, pointLat));
        }
      }
    }
  }



  @Override
  public HealpixNestedFixedRadiusConeComputer newComputer() {
    return new NestedLargeCell(this.hn, this.rRad);
  }

  private static double normalizeLon(final double lonRad) {
    if (lonRad < 0 || TWO_PI < lonRad) {
      final int k = HackersDelight.floorInt(lonRad * ONE_OVER_2PI);
      return lonRad - k * TWO_PI;
    }
    return lonRad;
  }

  private static double computeDeltaLon(final double lonARad, final double lonBRad) {
    final double res = abs(lonARad - lonBRad);
    return res <= PI ? res : TWO_PI - res;
  }

  private final double squareOfSinOfHalfR(final double deltaLonRad, final double deltaLatRad,
      final double cosLat, final double cosLatCenter) {
    return angDistComputer.squareOfsinOfhalfDistInRad2(deltaLonRad, deltaLatRad, cosLat, cosLatCenter);
  }

  private final double squareOfSinOfHalfOf(final double angleRad) {
    final double sin = angDistComputer.sin(0.5 * angleRad);
    return sin * sin;
  }
  
  private final double squareOfSinOfHalfOfDeltaLon(final double deltaLonRad) {
    final double sin = Math.sin(0.5 * deltaLonRad);
    return sin * sin;
  }

  /**
   * Return a value between 0 and 1, or NaN if deltaLatRad is larger than the radius of the cone.
   * We recall that any comparison (<, >, <=, >=, ==) other than (!=)  with NaN returns false 
   * @param deltaLatRad
   * @param cosLat
   * @param cosLatCenter
   * @return
   */
  private final double squareOfSinOfHalfDeltaLon(final double deltaLatRad,
      final double cosLat, final double cosLatCenter) {
    assert cosLat >= 0 && cosLatCenter >= 0;
    final double n = this.squareOfsinOfHalfR - squareOfSinOfHalfOf(deltaLatRad);
    final double d = cosLat * cosLatCenter;
    if (n < 0) { //  || !isFinite(n) || !isFinite(d)
      return Double.NaN;
    } else if (n >= d) {
      return 1; 
    } else {
      return n / d;  // return NaN if n or d is NaN
    }
  }

}
