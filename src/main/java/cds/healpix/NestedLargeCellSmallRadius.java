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

import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.CompassPoint.MainWind;

import cds.healpix.common.math.HackersDelight;

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
import static cds.healpix.common.math.Math.isFinite;

final class NestedLargeCellSmallRadius implements HealpixNestedFixedRadiusConeComputer {

  // private static final Logger LOG = LoggerFactory.getLogger("fr.unistra.cds.healpix.cone"); // NestedLargeCell.class.getName()

  private final int nIterMax;
  private final double epsLat;

  private final HealpixNested hn;
  private final double rRad;
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

  public NestedLargeCellSmallRadius(final HealpixNested hnTarget, final double radiusRad) {
    this(hnTarget, radiusRad, radiusRad * 1e-4, 10); 
  }

  protected NestedLargeCellSmallRadius(final HealpixNested hnTarget, final double radiusRad,
      final double epsLat, final int nIterMax) {
    assert radiusRad > 0;
    this.hn = hnTarget;
    this.epsLat = epsLat;
    this.nIterMax = nIterMax;
    this.rRad = radiusRad;
    this.angDistComputer = AngularDistanceComputer.getComputer(this.rRad);
    final double sinHalfR = angDistComputer.sin(0.5 * radiusRad); 
    //    this.twoSineOfHalfConeRadius = 2 * sinHalfR;
    this.squareOfsinOfHalfR = sinHalfR * sinHalfR;

    this.hComputer = this.hn.newHashComputer();
    this.neigSelector = this.hn.newNeighbourSelector();
    this.vpComputer = this.hn.newVerticesAndPathComputer();
    this.neigList = new NeighbourList(-1);
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
    assert iNeig == this.neigList.size() + 1;
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

    coneCenterLonRad = normalizeLon(coneCenterLonRad);
    assert 0 <= coneCenterLonRad && coneCenterLonRad <= TWO_PI;

    // Compute hash of the cell containing the cone center
    final long centerHash = this.hComputer.hash(coneCenterLonRad, coneCenterLatRad); // lat is checked here!
    assert -HALF_PI <= coneCenterLatRad && coneCenterLatRad <= HALF_PI;

    // Get the coordinates of the four vertices
    this.vpComputer.vertices(centerHash, this.vertices);

    final double latMax = coneCenterLatRad + this.rRad;
    final double latMin = coneCenterLatRad - this.rRad;

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
      if (squareOfSinOfHalfOf(computeDeltaLon(coneCenterLonRad, this.vertices.get(E)[LON_INDEX]))
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
      if (squareOfSinOfHalfOf(computeDeltaLon(coneCenterLonRad, this.vertices.get(W)[LON_INDEX]))
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


  // Use locally the orthographic (sinus) projection
  private double toLocalX(final double lonRad, final double latRad,
      final double coneCenterLonRad, final double coneCenterLatRad, final double cosConeCenterLat) {
    // return smallCos(latRad - coneCenterLatRad) * smallSin(lonRad - coneCenterLonRad);
    return cosConeCenterLat * smallSin(lonRad - coneCenterLonRad) / smallCos(lonRad - coneCenterLonRad);
  }

  private double toLocalY(final double latRad, final double coneCenterLatRad) {
    return smallSin(latRad - coneCenterLatRad);
  }

  private double straitLineSlope(double xA, double yA, double xB, double yB) {
    return (yB - yA) / (xB - xA); 
  }

  private double straitLineIntercept(double xA, double yA, double xB, double yB) {
    return (yA  * xB - xA * yB) / (xB - xA); 
  }

  private double smallSin(final double aRad) {
    return this.angDistComputer.sin(aRad);
    // assert aRad < 1e-5; // ~ 2arcsec
    // return aRad - ONE_OVER_6 * aRad * aRad * aRad; // + O(x^-25)
    // For a better precision, we do not use x * (1 + x * x / 6)
  }

  private double smallCos(final double aRad) {
    return this.angDistComputer.cos(aRad);
    // assert aRad < 1e-5; // ~ 2arcsec
    // final double x2 = aRad * aRad;
    // return 1 - 0.5 * x2 + ONE_OVER_24 * x2 * x2; // + O(x^20)
  }

  private double smallAsin(final double a) {
    return this.angDistComputer.asin(a);
  }

  private double latOfParallelLineTangentToCone(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double lonA, final double latA, final double lonB, final double latB,
      final boolean northSolution,
      final double discriminantToleranceInRadiusUnits) {
    // First project using the orthographic (SIN) projection
    // we do not use the ARC projection to avoid computing angular distances and cos.
    // With sin, we can work only on DeltaLon, DeltaLat, which allow faster trigonometric
    // computation (for small angles).
    final double xA = toLocalX(lonA, latA, coneCenterLonRad, coneCenterLatRad, cosConeCenterLat);
    final double yA = toLocalY(latA, coneCenterLatRad);
    final double xB = toLocalX(lonB, latB, coneCenterLonRad, coneCenterLatRad, cosConeCenterLat);
    final double yB = toLocalY(latB, coneCenterLatRad); 
    final double slope = straitLineSlope(xA, yA, xB, yB);
    final double r = smallSin(this.rRad); // The center of the circle is (0, 0)
    // LOG.debug("xA: {}; yA: {}; xB: {}; yB: {}; slope: {}; r: {}", xA, yA, xB, yB, slope, r);
    if (!isFinite(slope)) { // Vertical line
      assert abs(xB - xA) < 1e-15; // xA == xB (numerically)
      if (abs(xA) < (r * (1 + discriminantToleranceInRadiusUnits))) {
        return coneCenterLatRad;
      } else {
        return Double.NaN;
      }
    } else { // Normal line
      double intercept = straitLineIntercept(xA, yA, xB, yB);
      double a = slope * slope + 1; // x^2 coef. of the quadratic eq. of cone/strait line intersection
      double b = 2 * slope * intercept; // x coef. of the quadratic eq. of cone/strait line intersection
      double c = intercept  * intercept - r * r; // cte coef. of the quadratic eq. of cone/strait line intersection
      boolean hasSolution = (b * b - 4 * a * c) >= -r * discriminantToleranceInRadiusUnits;
      // LOG.debug("intercept: {}; a: {}; b: {}; c: {}; discri: {}", intercept, a, b, c,  (b * b - 4 * a * c));
      if (hasSolution) {
        double interceptUniqSolution = (northSolution ? r : -r) * Math.sqrt(a);
        double bUniqSolution = 2 * slope * interceptUniqSolution;
        double x = -bUniqSolution / (2 * a);
        double y = slope * x + interceptUniqSolution;
        // orthographic (SIN) deprojection:
        //   x3d = sqrt(1 - x * x - y * y);
        //   y3d = x;
        //   z3d = y;
        // LOG.debug("iUniqSol: {}; x: {}; y: {};", interceptUniqSolution, x, y);
        return coneCenterLatRad + smallAsin(y);
      } else { // NE edge do no intercept the cone (in the flat sky approximation)
        return Double.NaN;
      }
    }
  }

  private final void edgesNEWinEQR(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testNE, final boolean testNW) {    
    final double[] lonLatN = this.vertices.get(N);
    final double[] lonLatE = this.vertices.get(E);
    // First solution in the flat sky approximation (can save loops with trigo functions)
    double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat,
        lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, 1);
    // LOG.debug("Test N-E/W edges in EQR. at: {} deg;", Math.toDegrees(latRad));
    if (!isFinite(latRad)) { return; }
    // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = NewtonMethod.newtonSolveEquatorialZone(latRad,  coneCenterLatRad, cosConeCenterLat,
        squareOfsinOfHalfR,  false, this.epsLat, this.nIterMax, this.angDistComputer);
    if (!isFinite(latRad)) { return; }
    final double zDxDyEq1 = sin(latRad);
    final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        latRad - coneCenterLatRad, zDxDyEq1);
    // LOG.debug("Test N-E/W edges in EQR. zDxDyEq1: {}; lat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
    if (testNE) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, latRad));
    }
    if (testNW) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, latRad));
    } 
  }

  private final void edgesSEWinEQR(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testSE, final boolean testSW) {
    final double[] lonLatS = this.vertices.get(S);
    final double[] lonLatE = this.vertices.get(E);
    // First solution in the flat sky approximation (can save loops with trigo functions)
    double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat, 
        lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], false, 1);
    // LOG.debug("Test S-E/W edges in EQR. at: {} deg;", Math.toDegrees(latRad));
    if (!isFinite(latRad)) { return; }
    // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = NewtonMethod.newtonSolveEquatorialZone(latRad,  coneCenterLatRad, cosConeCenterLat,
        squareOfsinOfHalfR,  true, this.epsLat, this.nIterMax, this.angDistComputer);
    if (!isFinite(latRad)) { return; }
    final double zDxDyEq1 = sin(latRad);
    final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    // LOG.debug("Test S-E/W edge in EQR. zDxDyEq1: {}; lat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
    if (testSE) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, latRad));
    }
    if (testSW) {
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, latRad));
    }
  }

  private final void edgesNEWinNPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testNE, final boolean testNW, final boolean testN) {
    // Remark: North pole not included, else (isNotSelectedNE || isNotSelectedNW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double[] lonLatN = this.vertices.get(N);
    final double[] lonLatE = this.vertices.get(E);
    final double[] lonLatW = this.vertices.get(W);
    if (testNE || testN) {
      // First solution in the flat sky approximation (can save loops with trigo functions)
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, 1);
      // LOG.debug("Test NE edge in NPC. lat: {} deg;", Math.toDegrees(latRad));
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, this.epsLat, this.nIterMax,
          this.angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      // LOG.debug("Test NE edge in NPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, latRad));
      if (coneCenterLonModHalfPi + deltaLon > PI_OVER_FOUR && lonLatN[LAT_INDEX] < HALF_PI) { // Change the base pixe (exept at the north pole)
        // It means that the NE edge is vertical
        // First solution in the flat sky approximation (can save loops with trigo functions)
        // Use the symmetric of the W vertex by the vertical NE axis
        latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat,
            lonLatN[LON_INDEX], lonLatN[LAT_INDEX], 2 * lonLatN[LON_INDEX] - lonLatW[LON_INDEX], lonLatW[LAT_INDEX],
            true, 1);
        // LOG.debug("Test NE edges in NPC bis. lat: {} deg;", Math.toDegrees(latRad));
        // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi - HALF_PI,
            coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, this.epsLat, this.nIterMax,
            this.angDistComputer);
        zDxDyEq1 = sin(latRad);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        // LOG.debug("Test NE edge in NPC bis. zDxDyEq1: {}; lat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
        if (isFinite(deltaLon)) {
          insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, latRad));
        }
      }
    }
    if (testNW || testN) {
      // First solution in the flat sky approximation (can save loops with trigo functions)
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatN[LON_INDEX], lonLatN[LAT_INDEX], true, 1);
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, this.epsLat, this.nIterMax,
          this.angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      // LOG.debug("Test NW edge in NPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, latRad));
      if (coneCenterLonModHalfPi - deltaLon < 0 && lonLatN[LAT_INDEX] < HALF_PI) { // Change the base pixel!
        // First solution in the flat sky approximation (can save loops with trigo functions)
        // Use the symmetric of the E vertex by the vertical NW axis
        latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat,
            2 * lonLatN[LON_INDEX] - lonLatE[LON_INDEX], lonLatE[LAT_INDEX], lonLatN[LON_INDEX], lonLatN[LAT_INDEX],
            true, 1);
        // LOG.debug("Test NW edges in NPC bis. lat: {} deg;", Math.toDegrees(latRad));
        // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, coneCenterLonModHalfPi + HALF_PI,
            coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, this.epsLat, this.nIterMax,
            this.angDistComputer);
        zDxDyEq1 = sin(latRad);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        // E.g. in Aladin: draw circle(089.91417, +47.04352, 59')
        // LOG.debug("Test NW edge in NPC bis. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
        if (isFinite(deltaLon)) {
          insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, latRad));
        }
      }
    }
  }

  private final void edgesNEWinSPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testNE, final boolean testNW) {
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double[] lonLatN = this.vertices.get(N);
    final double[] lonLatE = this.vertices.get(E);
    final double[] lonLatW = this.vertices.get(W);
    // Remark: North pole not included, else (isNotSelectedNE || isNotSelectedNW) == false
    if (testNE) {
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, 1);
      if (!isFinite(latRad)) { return; }
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, this.epsLat, this.nIterMax,
          this.angDistComputer);
      if (!isFinite(latRad)) { return; }
        final double zDxDyEq1 = sin(latRad);
        final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        // LOG.debug("Test NE edge in SPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, latRad));
    }
    if (testNW) {
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatN[LON_INDEX], lonLatN[LAT_INDEX], true, 1);
      if (!isFinite(latRad)) { return; }
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, this.epsLat, this.nIterMax,
          this.angDistComputer);
      if (!isFinite(latRad)) { return; }
      final double zDxDyEq1 = sin(latRad);
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      // LOG.debug("Test NW edge in SPC. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, latRad));
    }
  }

  private final void edgesSEWinNPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testSE, final boolean testSW) {
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double[] lonLatS = this.vertices.get(S);
    final double[] lonLatE = this.vertices.get(E);
    final double[] lonLatW = this.vertices.get(W);
    // Remark: South pole not included, else (testSE || testSW) == false
    if (testSE) {
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], false, 1);
      if (!isFinite(latRad)) { return; }
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, this.epsLat, this.nIterMax,
          this.angDistComputer);
      if (!isFinite(latRad)) { return; }
      final double zDxDyEq1 = sin(latRad);
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      // LOG.debug("Test SE edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, latRad));
    }
    if (testSW) {
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatS[LON_INDEX], lonLatS[LAT_INDEX], false, 1);
      if (!isFinite(latRad)) { return; }
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, this.epsLat, this.nIterMax,
          this.angDistComputer);
      if (!isFinite(latRad)) { return; }
      final double zDxDyEq1 = sin(latRad);
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      // LOG.debug("Test SW edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
      insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, latRad));
    }
  }

  private final void edgesSEWinSPC(
      final double coneCenterLonRad, final double coneCenterLatRad,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final boolean testSE, final boolean testSW, final boolean testS) {
    // Remark: South pole not included, else (testSE || testSW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double[] lonLatS = this.vertices.get(S);
    final double[] lonLatE = this.vertices.get(E);
    final double[] lonLatW = this.vertices.get(W);
    if (testSE || testS) {
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, 1); // A REVOIR!!!
      // LOG.debug("Test S-E/W edges in SPC. at: {} deg;", Math.toDegrees(latRad));
      // assert latRad < 0;
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, this.epsLat, this.nIterMax,
          this.angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      if (coneCenterLonModHalfPi + deltaLon > PI_OVER_FOUR && lonLatS[LAT_INDEX] > -HALF_PI) {
        latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat,
            lonLatS[LON_INDEX], lonLatS[LAT_INDEX], 2 * lonLatS[LON_INDEX] - lonLatW[LON_INDEX], lonLatW[LAT_INDEX],
            false, 1);
        // LOG.debug("Test S-E/W edges in SPC v2. at: {} deg;", Math.toDegrees(latRad));
        // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi - HALF_PI,
            -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, this.epsLat, this.nIterMax,
            this.angDistComputer);
        zDxDyEq1 = sin(latRad);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            coneCenterLatRad - latRad, zDxDyEq1);
      }
      // LOG.debug("Test SE edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
      if (isFinite(latRad)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad + deltaLon, latRad));
      }
    }
    if (testSW || testS) {
      double latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatS[LON_INDEX], lonLatS[LAT_INDEX], false, 1);
      // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, this.epsLat, this.nIterMax,
          this.angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      if (coneCenterLonModHalfPi - deltaLon < 0 && lonLatS[LAT_INDEX] > -HALF_PI) {
        latRad = latOfParallelLineTangentToCone(coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat, 
            2 * lonLatS[LON_INDEX] - lonLatE[LON_INDEX], lonLatE[LAT_INDEX], lonLatS[LON_INDEX], lonLatS[LAT_INDEX],
            false, 1);
        // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, coneCenterLonModHalfPi + HALF_PI,
            -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, this.epsLat, this.nIterMax,
            this.angDistComputer);
        zDxDyEq1 = sin(latRad);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            coneCenterLatRad - latRad, zDxDyEq1);
      }
      // LOG.debug("Test SW edge in SPC/EQR. zDxDyEq1: {}; pointLat: {} deg; deltaLon: {} deg", zDxDyEq1, Math.toDegrees(latRad), Math.toDegrees(deltaLon));
      if (isFinite(latRad)) {
        insertSortRmDuplicates(this.hComputer.hash(coneCenterLonRad - deltaLon, latRad));
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
