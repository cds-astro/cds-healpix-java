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

import static cds.healpix.CompassPoint.Cardinal.E;
import static cds.healpix.CompassPoint.Cardinal.N;
import static cds.healpix.CompassPoint.Cardinal.S;
import static cds.healpix.CompassPoint.Cardinal.W;
import static cds.healpix.Healpix.isLatInNorthPolarCap;
import static cds.healpix.Healpix.isLatInSouthPolarCap;
import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.abs;
import static cds.healpix.common.math.Math.asin;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.isFinite;

import java.util.EnumMap;

import cds.healpix.CompassPoint.Cardinal;

final class SmallConeOrdinalHashComputer implements ConeOrdinalHashComputer {

  public static final SmallConeOrdinalHashComputer UI = new SmallConeOrdinalHashComputer();
  
  private SmallConeOrdinalHashComputer() { }
  
  @Override
  public int computeOrdinalHash(final double coneCenterLonRad, final double coneCenterLatRad,
      final double coneRadiusRad, final HashComputer hashComputer, final AngularDistanceComputer angDistComputer, 
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities
      final long hashCenterAtSmallestDepth, final VerticesAndPathComputer verticesComputerAtSmallestDepth,
      final EnumMap<Cardinal, double[]> vertices,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      // Result params
      final long[] result) {
    final double sinOfConeRadius = angDistComputer.sin(coneRadiusRad);
    verticesComputerAtSmallestDepth.vertices(hashCenterAtSmallestDepth, vertices);
    int resultSize = 0;
    // Compute N-E/W hash
    if (isLatInNorthPolarCap(coneCenterLatRad)) {
      resultSize = edgesNEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else if (isLatInSouthPolarCap(coneCenterLatRad + coneRadiusRad)) {
      resultSize = edgesNEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else {
      resultSize = edgesNEWinEQR(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    }
    // Compute S-E/W hash
    if (isLatInNorthPolarCap(coneCenterLatRad - coneRadiusRad)) {   // SE/SW edges in north polar cap
      resultSize = edgesSEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else if (isLatInSouthPolarCap(coneCenterLatRad)) { // SE/SW edges in south polar cap
      resultSize = edgesSEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else { // SE/SW edges in the EQR
      resultSize = edgesSEWinEQR(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    }
    return resultSize;
  }

 //Use locally the orthographic (sinus) projection
 private double toLocalX(final double lonRad, final double latRad,
     final double coneCenterLonRad, final double coneCenterLatRad, final double cosConeCenterLat,
     final AngularDistanceComputer angDistComputer) {
   return cosConeCenterLat * angDistComputer.sin(lonRad - coneCenterLonRad)
       / angDistComputer.cos(lonRad - coneCenterLonRad);
 }

 private static double toLocalY(final double latRad, final double coneCenterLatRad,
     final AngularDistanceComputer angDistComputer) {
   return angDistComputer.sin(latRad - coneCenterLatRad);
 }

 private double latOfParallelLineTangentToCone(
     final double sinOfConeRadius,
     final double coneCenterLonRad, final double coneCenterLatRad,
     final double cosConeCenterLat, final double sinConeCenterLat,
     final double lonA, final double latA, final double lonB, final double latB,
     final boolean northSolution,
     final AngularDistanceComputer angDistComputer) {
   // First project using the orthographic (SIN) projection
   // we do not use the ARC projection to avoid computing angular distances and cos.
   // With sin, we can work only on DeltaLon, DeltaLat, which allow faster trigonometric
   // computation for small angles.
   final double xA = toLocalX(lonA, latA, coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, angDistComputer);
   final double yA = toLocalY(latA, coneCenterLatRad, angDistComputer);
   final double xB = toLocalX(lonB, latB, coneCenterLonRad, coneCenterLatRad, cosConeCenterLat, angDistComputer);
   final double yB = toLocalY(latB, coneCenterLatRad, angDistComputer); 
   final double slope = (yB - yA) / (xB - xA);
   final double r = sinOfConeRadius; // The center of the circle is (0, 0)
// LOG.debug("xA: {}; yA: {}; xB: {}; yB: {}; slope: {}; r: {}", xA, yA, xB, yB, slope, r);
   if (!isFinite(slope)) { // Vertical line
     assert abs(xB - xA) < 1e-15; // xA == xB (numerically)
     return coneCenterLatRad;
   } else { // Normal line
     final double a = slope * slope + 1; // x^2 coef. of the quadratic eq. of cone/strait line intersection
     final double interceptUniqSolution = (northSolution ? r : -r) * Math.sqrt(a);
     final double bUniqSolution = 2 * slope * interceptUniqSolution;
     final double x = -bUniqSolution / (2 * a);
     final double y = slope * x + interceptUniqSolution;
       // orthographic (SIN) deprojection:
       //   x3d = sqrt(1 - x * x - y * y);
       //   y3d = x;
       //   z3d = y;
// LOG.debug("iUniqSol: {}; x: {}; y: {};", interceptUniqSolution, x, y);
       return coneCenterLatRad + angDistComputer.asin(y);
   }
 }
  
  private final int edgesNEWinEQR(
      final double coneCenterLonRad, final double coneCenterLatRad, final double coneRadiusRad,
      final HashComputer hashComputer, final AngularDistanceComputer angDistComputer,
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities 
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      final double sinOfConeRadius, final EnumMap<Cardinal, double[]> vertices,
      // Result params
      final long[] result, int resultSize) {
    final double[] lonLatN = vertices.get(N);
    final double[] lonLatE = vertices.get(E);
    double latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat,
        lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, angDistComputer);
// LOG.debug("Test N-E/W edges in EQR. at: {} deg;", Math.toDegrees(latRad));
    // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = NewtonMethod.newtonSolveEquatorialZone(latRad,  coneCenterLatRad, cosConeCenterLat,
        squareOfsinOfHalfR,  false, relativePrecision, nIterMax, angDistComputer);
    final double zDxDyEq1 = sin(latRad);
    if (isLatInNorthPolarCap(latRad)) {
      resultSize = edgesNEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else if (isLatInSouthPolarCap(latRad)) {
      resultSize = edgesNEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else {
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ?  hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
      result[resultSize++] = isFinite(deltaLon) ?  hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1L;
    }
    return resultSize;
  }

  private final int edgesSEWinEQR(
      final double coneCenterLonRad, final double coneCenterLatRad, final double coneRadiusRad,
      final HashComputer hashComputer, final AngularDistanceComputer angDistComputer,
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      final double sinOfConeRadius, final EnumMap<Cardinal, double[]> vertices,
      // Result params
      final long[] result, int resultSize) {
    final double[] lonLatS = vertices.get(S);
    final double[] lonLatE = vertices.get(E);
    // First solution in the flat sky approximation (can save loops with trigo functions)
    double latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat, 
        lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], false, angDistComputer);
// LOG.debug("Test S-E/W edges in EQR. at: {} deg;", Math.toDegrees(latRad));
    // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = NewtonMethod.newtonSolveEquatorialZone(latRad,  coneCenterLatRad, cosConeCenterLat,
        squareOfsinOfHalfR,  true, relativePrecision, nIterMax, angDistComputer);
    final double zDxDyEq1 = sin(latRad);
    if (isLatInNorthPolarCap(latRad)) {
      resultSize = edgesSEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else if (isLatInSouthPolarCap(latRad)) {
      resultSize = edgesSEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          sinOfConeRadius, vertices, result, resultSize);
    } else {
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1L;
    }
    return resultSize;
  }

  private final int edgesNEWinNPC(
      final double coneCenterLonRad, final double coneCenterLatRad, final double coneRadiusRad,
      final HashComputer hashComputer, final AngularDistanceComputer angDistComputer,
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      final double sinOfConeRadius, final EnumMap<Cardinal, double[]> vertices,
      // Result params
      final long[] result, int resultSize) {
    // Remark: North pole not included, else (isNotSelectedNE || isNotSelectedNW) == false
    final double[] lonLatN = vertices.get(N);
    final double[] lonLatE = vertices.get(E);
    final double[] lonLatW = vertices.get(W);
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final boolean containsNorthPole = coneCenterLatRad + coneRadiusRad > HALF_PI;
    if (containsNorthPole) {
      // Hypothese, the opposite cell (other side of the pole) never contains a perculier point
      final double dCN = HALF_PI - coneCenterLatRad;
      final double d1 = coneRadiusRad + dCN;
      final double d2 = coneRadiusRad - dCN;
      final double cosL = Math.cos(coneCenterLonModHalfPi);
      final double sinL = Math.sin(coneCenterLonModHalfPi);
      // Look at the south east value in the west celllatRad = asin(zDxDyEq1);
      double latStart = 0.5 * (d1 * cosL + d2 * sinL);
      // Exact calculation.  TODO: use firs the flat sky approximations (with y>1, so need to look at x, ...
      double latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latStart, true, coneCenterLonModHalfPi - HALF_PI,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
          angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1;
      // Look at the south west value in the west cell
      latStart =  0.5 * (d2 * cosL + d1 * sinL);
      // Exact calculation.  TODO: use firs the flat sky approximations (with y>1, so need to look at x, ...
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latStart, false, HALF_PI + coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ?  hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1;
    } else {
      // North East
      // - First solution in the flat sky approximation (can save loops with trigo functions)
      double latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, angDistComputer);
      // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
          angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
      } else { // Peculiar point in the NE neighbour base cell
        // - First solution in the flat sky approximation (can save loops with trigo functions)
        latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat,
            lonLatN[LON_INDEX], lonLatN[LAT_INDEX], 2 * lonLatN[LON_INDEX] - lonLatW[LON_INDEX], lonLatW[LAT_INDEX],
            true, angDistComputer);
        // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi - HALF_PI,
            coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
            angDistComputer);
        zDxDyEq1 = sin(latRad);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the North pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
        } else { // Near from the North pole => compute south east
          //  - First solution in the flat sky approximation (can save loops with trigo functions)
          latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
              cosConeCenterLat, sinConeCenterLat,
              lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], false, angDistComputer);
          // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
          latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi - HALF_PI,
              coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
              angDistComputer);
          zDxDyEq1 = sin(latRad);
          deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
              coneCenterLatRad - latRad, zDxDyEq1);
          result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
        }
      }
      // North West
      // - First solution in the flat sky approximation (can save loops with trigo functions)
      latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatN[LON_INDEX], lonLatN[LAT_INDEX], true, angDistComputer);
      // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
      } else { // Peculiar point in the NW neighbour base cell
        // - First solution in the flat sky approximation (can save loops with trigo functions)
        latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat,
            2 * lonLatN[LON_INDEX] - lonLatE[LON_INDEX], lonLatE[LAT_INDEX], lonLatN[LON_INDEX], lonLatN[LAT_INDEX],
            true, angDistComputer);
        // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, HALF_PI + coneCenterLonModHalfPi,
            coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
            angDistComputer);
        zDxDyEq1 = sin(latRad);
        latRad = asin(zDxDyEq1);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the North pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
        } else { // Near from the North pole => compute south west
          // - First solution in the flat sky approximation (can save loops with trigo functions)
          latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
              cosConeCenterLat, sinConeCenterLat,
              lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatN[LON_INDEX], lonLatN[LAT_INDEX], false, angDistComputer);
          // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
          latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, HALF_PI + coneCenterLonModHalfPi,
              coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
              angDistComputer);
          zDxDyEq1 = sin(latRad);
          deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
              coneCenterLatRad - latRad, zDxDyEq1);
          result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1L;
        }
      }
    }
    return resultSize;
  }

  private final int edgesSEWinNPC(
      final double coneCenterLonRad, final double coneCenterLatRad, final double coneRadiusRad,
      final HashComputer hashComputer, final AngularDistanceComputer angDistComputer,
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      final double sinOfConeRadius, final EnumMap<Cardinal, double[]> vertices,
      // Result params
      final long[] result, int resultSize) {
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double[] lonLatS = vertices.get(S);
    final double[] lonLatE = vertices.get(E);
    final double[] lonLatW = vertices.get(W);
    // South-east
    double latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat,
        lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], false, angDistComputer);
    // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi,
        coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
        angDistComputer);
    double zDxDyEq1 = sin(latRad);
    double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
    } else { // Look at the neighbour base cell, compute south east (vertical edge, i.e. // to NE)
      final double[] lonLatN = vertices.get(S);
      latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, angDistComputer);
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, true, coneCenterLonModHalfPi - HALF_PI,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
    }
    // South-west
    latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat,
        lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatS[LON_INDEX], lonLatS[LAT_INDEX], false, angDistComputer);
    // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, coneCenterLonModHalfPi,
        coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
        angDistComputer);
    zDxDyEq1 = sin(latRad);
    deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
    } else { // Look at the neighbour base cell, comute south west (vertical edge, i.e. // NW)
      final double[] lonLatN = vertices.get(S);
      latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatW[LON_INDEX], lonLatW[LAT_INDEX], true, angDistComputer);
      latRad = NewtonMethod.newtonSolveNorthPolarCapZone(latRad, false, HALF_PI + coneCenterLonModHalfPi,
          coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1L;
    }
    return resultSize;
  }



  private final int edgesSEWinSPC(final double coneCenterLonRad, final double coneCenterLatRad, final double coneRadiusRad,
      final HashComputer hashComputer, final AngularDistanceComputer angDistComputer,
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      final double sinOfConeRadius, final EnumMap<Cardinal, double[]> vertices,
      // Result params
      final long[] result, int resultSize) {
    final double[] lonLatS = vertices.get(S);
    final double[] lonLatE = vertices.get(E);
    final double[] lonLatW = vertices.get(W);
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final boolean containsSouthPole = coneRadiusRad - cosConeCenterLat > HALF_PI;
    if (containsSouthPole) {
      final double dCS = HALF_PI + coneCenterLatRad;
      final double d1 = coneRadiusRad + dCS;
      final double d2 = coneRadiusRad - dCS;
      final double cosL = Math.cos(coneCenterLonModHalfPi);
      final double sinL = Math.sin(coneCenterLonModHalfPi);
      // Look at the south east value in the west celllatRad = asin(zDxDyEq1);
      double latStart = 0.5 * (d1 * cosL + d2 * sinL);
      // Exact calculation.  TODO: use firs the flat sky approximations (with y>1, so need to look at x, ...
      double latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(latStart, true, coneCenterLonModHalfPi - HALF_PI,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
          angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1;
      // Look at the south west value in the west cell
      latStart =  0.5 * (d2 * cosL + d1 * sinL);
      // Exact calculation.  TODO: use firs the flat sky approximations (with y>1, so need to look at x, ...
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(latStart, false, HALF_PI + coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1;
    } else {
      // South east
      // - First solution in the flat sky approximation (can save loops with trigo functions)
      double latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], false, angDistComputer);
      // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
          angDistComputer);
      double zDxDyEq1 = sin(latRad);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
      } else { // Peculiar point in the SE neighbour base cell
        // - First solution in the flat sky approximation (can save loops with trigo functions)
        latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat,
            lonLatS[LON_INDEX], lonLatS[LAT_INDEX], 2 * lonLatS[LON_INDEX] - lonLatW[LON_INDEX], lonLatW[LAT_INDEX],
            false, angDistComputer);
        // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi - HALF_PI,
            -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
            angDistComputer);
        zDxDyEq1 = sin(latRad);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            coneCenterLatRad - latRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the South pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
        } else { // Near from the South pole => compute north east
          // - First solution in the flat sky approximation (can save loops with trigo functions)
          latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
              cosConeCenterLat, sinConeCenterLat,
              lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, angDistComputer);
          // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
          latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi - HALF_PI,
              -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
              angDistComputer);
          zDxDyEq1 = sin(latRad);
          deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
              coneCenterLatRad - latRad, zDxDyEq1);
          result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
        }
      }
      // South west
      // - First solution in the flat sky approximation (can save loops with trigo functions)
      latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatS[LON_INDEX], lonLatS[LAT_INDEX], false, angDistComputer);
      // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
      } else { // Peculiar point in the WE neighbour base cell
        // - First solution in the flat sky approximation (can save loops with trigo functions)
        latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
            cosConeCenterLat, sinConeCenterLat,
            2 * lonLatS[LON_INDEX] - lonLatE[LON_INDEX], lonLatE[LAT_INDEX], lonLatS[LON_INDEX], lonLatS[LAT_INDEX],
            false, angDistComputer);
        // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
        latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, HALF_PI + coneCenterLonModHalfPi,
            -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
            angDistComputer);
        zDxDyEq1 = sin(latRad);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the North pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
        } else { // Near from the North pole => compute south west
          // - First solution in the flat sky approximation (can save loops with trigo functions)
          latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
              cosConeCenterLat, sinConeCenterLat,
              lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatS[LON_INDEX], lonLatS[LAT_INDEX], true, angDistComputer);
          // - Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
          latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, HALF_PI + coneCenterLonModHalfPi,
              -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
              angDistComputer);
          zDxDyEq1 = sin(latRad);
          deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
              coneCenterLatRad - latRad, zDxDyEq1);
          result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1L;
        }
      }
    }
    return resultSize;
  }

  private final int edgesNEWinSPC(final double coneCenterLonRad, final double coneCenterLatRad, final double coneRadiusRad,
      final HashComputer hashComputer, final AngularDistanceComputer angDistComputer,
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      final double sinOfConeRadius, final EnumMap<Cardinal, double[]> vertices,
      // Result params
      final long[] result, int resultSize) {
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double[] lonLatN = vertices.get(N);
    final double[] lonLatE = vertices.get(E);
    final double[] lonLatW = vertices.get(W);
    // North-east
    double latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat,
        lonLatN[LON_INDEX], lonLatN[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], true, angDistComputer);
    // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi,
        -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
        angDistComputer);
    double zDxDyEq1 = sin(latRad);
    double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
    } else { // Look at the neighbour base cell, comute north east
      final double[] lonLatS = vertices.get(S);
      latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatE[LON_INDEX], lonLatE[LAT_INDEX], false, angDistComputer);
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, true, coneCenterLonModHalfPi - HALF_PI,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, false, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
    }
    // North-west
    latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
        cosConeCenterLat, sinConeCenterLat,
        lonLatW[LON_INDEX], lonLatW[LAT_INDEX], lonLatN[LON_INDEX], lonLatN[LAT_INDEX], true, angDistComputer);
    // Then exact calculation (hopefully only 1 iteration if the flat sky approx is good)
    latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, coneCenterLonModHalfPi,
        -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
        angDistComputer);
    zDxDyEq1 = sin(latRad);
    deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
    } else { // Look at the neighbour base cell, comute south west
      final double[] lonLatS = vertices.get(S);
      latRad = latOfParallelLineTangentToCone(sinOfConeRadius, coneCenterLonRad, coneCenterLatRad,
          cosConeCenterLat, sinConeCenterLat,
          lonLatS[LON_INDEX], lonLatS[LAT_INDEX], lonLatW[LON_INDEX], lonLatW[LAT_INDEX], false, angDistComputer);
      latRad = -NewtonMethod.newtonSolveNorthPolarCapZone(-latRad, false, HALF_PI + coneCenterLonModHalfPi,
          -coneCenterLatRad, cosConeCenterLat, squareOfsinOfHalfR, true, relativePrecision, nIterMax,
          angDistComputer);
      zDxDyEq1 = sin(latRad);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1L;
    }
    return resultSize;
  }

}
