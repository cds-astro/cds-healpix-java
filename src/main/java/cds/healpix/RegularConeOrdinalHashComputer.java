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

import static cds.healpix.Healpix.isLatInNorthPolarCap;
import static cds.healpix.Healpix.isLatInSouthPolarCap;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.asin;
import static cds.healpix.common.math.Math.isFinite;

import java.util.EnumMap;

import cds.healpix.CompassPoint.Cardinal;

final class RegularConeOrdinalHashComputer implements ConeOrdinalHashComputer {

  public static final RegularConeOrdinalHashComputer UI = new RegularConeOrdinalHashComputer();
  
  private RegularConeOrdinalHashComputer() { }
  
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
    int resultSize = 0;
    // Compute N-E/W hash
    if (isLatInNorthPolarCap(coneCenterLatRad)) {
      resultSize = edgesNEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    } else if (isLatInSouthPolarCap(coneCenterLatRad + coneRadiusRad)) {
      resultSize = edgesNEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    } else {
      resultSize = edgesNEWinEQR(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    }
    // Compute S-E/W hash
    if (isLatInNorthPolarCap(coneCenterLatRad - coneRadiusRad)) {   // SE/SW edges in north polar cap
      resultSize = edgesSEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    } else if (isLatInSouthPolarCap(coneCenterLatRad)) { // SE/SW edges in south polar cap
      resultSize = edgesSEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    } else { // SE/SW edges in the EQR
      resultSize = edgesSEWinEQR(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    }
    return resultSize;
  }

  private final int edgesNEWinEQR(
      final double coneCenterLonRad, final double coneCenterLatRad, final double coneRadiusRad,
      final HashComputer hashComputer, final AngularDistanceComputer angDistComputer,
      // Algo parameters
      final double relativePrecision, final int nIterMax,
      // Pre-computed quantities
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      // Result params
      final long[] result, int resultSize) {
    final double zDxDyEq1 = NewtonMethod.newtonSolveEquatorialZone(
        Math.sin(coneCenterLatRad + coneRadiusRad * 0.9),
        sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
    final double latRad = asin(zDxDyEq1);
    if (isLatInNorthPolarCap(latRad)) {
      resultSize = edgesNEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    } else if (isLatInSouthPolarCap(latRad)) {
      resultSize = edgesNEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
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
      // Result params
      final long[] result, int resultSize) {
    final double zDxDyEq1 = NewtonMethod.newtonSolveEquatorialZone(
        Math.sin(coneCenterLatRad - coneRadiusRad * 0.9),
        sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
    final double latRad = asin(zDxDyEq1);
    if (isLatInNorthPolarCap(latRad)) {
      resultSize = edgesSEWinNPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    } else if (isLatInSouthPolarCap(latRad)) {
      resultSize = edgesSEWinSPC(coneCenterLonRad, coneCenterLatRad, coneRadiusRad,
          hashComputer, angDistComputer, relativePrecision, nIterMax, 
          cosConeCenterLat, sinConeCenterLat, twoSineOfHalfConeRadius, squareOfsinOfHalfR,
          result, resultSize);
    } else {
      final double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ?  hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
      result[resultSize++] = isFinite(deltaLon) ?  hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1L;
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
      // Result params
      final long[] result, int resultSize) {
    // Remark: North pole not included, else (isNotSelectedNE || isNotSelectedNW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final boolean containsNorthPole = coneCenterLatRad + coneRadiusRad > HALF_PI;
    if (containsNorthPole) {
      // Hypothese, the opposite cell (other side of the pole) never contains a perculier point
      final double dCN = HALF_PI - coneCenterLatRad;
      final double d1 = coneRadiusRad + dCN;
      final double d2 = coneRadiusRad - dCN;
      final double cosL = Math.cos(coneCenterLonModHalfPi);
      final double sinL = Math.sin(coneCenterLonModHalfPi);
      double latStart = 0.5 * (d1 * cosL + d2 * sinL);
      double zStart =  Math.sin(latStart);
      // look at the south east value in the east cell
      double zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, true, coneCenterLonModHalfPi - HALF_PI,
          sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
      double latRad = asin(zDxDyEq1);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1;
      // look at the south west value in the west cell
      latStart = 0.5 * (d2 * cosL + d1 * sinL);
      zStart =  Math.sin(latStart);
      zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, false, HALF_PI + coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ?  hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1;
    } else {
      // East
      final double zStart = Math.sin(coneCenterLatRad + coneRadiusRad * 0.9);
      double zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, true, coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
      // Remark: z can't be near from 1 since: distace min to pole = radius, and radius not small 
      double latRad = asin(zDxDyEq1);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
      } else { // Peculiar point in the NE neighbour base cell
        zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
            zStart, true, coneCenterLonModHalfPi - HALF_PI,
            sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
        latRad = asin(zDxDyEq1);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the North pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
        } else { // Near from the North pole => compute south east
          zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
              sinConeCenterLat, true, coneCenterLonModHalfPi - HALF_PI,
              sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
          latRad = asin(zDxDyEq1);
          deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
              coneCenterLatRad - latRad, zDxDyEq1);
          result[resultSize++] = isFinite(deltaLon) ?  hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
        }
      }
      // West
      zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, false, coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
      } else { // Peculiar point in the NW neighbour base cell
        zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
            zStart, false, HALF_PI + coneCenterLonModHalfPi,
            sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
        latRad = asin(zDxDyEq1);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the North pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
        } else { // Near from the North pole => compute south west 
          zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
              sinConeCenterLat, false, HALF_PI + coneCenterLonModHalfPi,
              sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
          latRad = asin(zDxDyEq1);
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
      // Result params
      final long[] result, int resultSize) {
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double zStart = Math.sin(coneCenterLatRad - relativePrecision * 0.9);
    // South-east
// System.out.println("zStart: " + zStart + "; coneCenterLonModHalfPi: " + coneCenterLonModHalfPi);
    double zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
        zStart, true, coneCenterLonModHalfPi,
        sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
    double latRad = asin(zDxDyEq1); 
    double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
// System.out.println("deltaLon: " + deltaLon + "; zDxDyEq1: " + zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
    } else { // Look at the neighbour base cell, comute south east
      zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          sinConeCenterLat, true, coneCenterLonModHalfPi - HALF_PI,
          sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
// System.out.println("deltaLon: " + deltaLon + "; zDxDyEq1: " + zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
    }
    // South-west
    zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
        zStart, false, coneCenterLonModHalfPi,
        sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
    latRad = asin(zDxDyEq1); 
    deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
    } else { // Look at the neighbour base cell, comute south west
      zDxDyEq1 = NewtonMethod.newtonSolveNorthPolarCapZone(
          sinConeCenterLat, false, HALF_PI + coneCenterLonModHalfPi,
          sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
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
      // Result params
      final long[] result, int resultSize) {
    // Remark: North pole not included, else (isNotSelectedNE || isNotSelectedNW) == false
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final boolean containsSouthPole = coneRadiusRad - cosConeCenterLat > HALF_PI;
    if (containsSouthPole) {
      // Hypothese, the opposite cell (other side of the pole) never contains a perculier point
      final double dSC = HALF_PI + coneCenterLatRad;
      final double d1 = coneRadiusRad + dSC;
      final double d2 = coneRadiusRad - dSC;
      final double cosL = Math.cos(coneCenterLonModHalfPi);
      final double sinL = Math.sin(coneCenterLonModHalfPi);
      double latStart = 0.5 * (d1 * cosL + d2 * sinL);
      double zStart =  Math.sin(latStart);
      // look at the south east value in the east cell
      double zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, true, coneCenterLonModHalfPi - HALF_PI,
          -sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
      double latRad = asin(zDxDyEq1);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1;
      // look at the south west value in the west cell
      latStart = 0.5 * (d2 * cosL + d1 * sinL);
      zStart =  Math.sin(latStart);
      zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, false, HALF_PI + coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad) : -1;
    } else {
      // South east
      final double zStart = Math.sin(-coneCenterLatRad + coneRadiusRad * 0.9);
      double zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, true, coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
      double latRad = asin(zDxDyEq1);
      double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
      } else { // Peculiar point in the SE neighbour base cell
        zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
            zStart, true, coneCenterLonModHalfPi - HALF_PI,
            -sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
        latRad = asin(zDxDyEq1);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            coneCenterLatRad - latRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the South pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
        } else { // Near from the South pole => compute north east
          zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
              -sinConeCenterLat, true, coneCenterLonModHalfPi - HALF_PI,
              -sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
          latRad = asin(zDxDyEq1);
          deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
              coneCenterLatRad - latRad, zDxDyEq1);
          result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : -1L;
        }
      }
      // South west
      zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          zStart, false, coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          latRad - coneCenterLatRad, zDxDyEq1);
      if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
        result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
      } else { // Peculiar point in the WE neighbour base cell
        zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
            zStart, false, HALF_PI + coneCenterLonModHalfPi,
            -sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
        latRad = asin(zDxDyEq1);
        deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
            latRad - coneCenterLatRad, zDxDyEq1);
        if (isFinite(deltaLon)) { // Far from the North pole
          result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
        } else { // Near from the North pole => compute south west 
          zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
              -sinConeCenterLat, false, HALF_PI + coneCenterLonModHalfPi,
              -sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
          latRad = asin(zDxDyEq1);
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
      // Result params
      final long[] result, int resultSize) {
    final double coneCenterLonModHalfPi = coneCenterLonRad % HALF_PI;
    final double zStart = Math.sin(-coneCenterLatRad + relativePrecision * 0.9);
    // North-east
    double zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
        zStart, true, coneCenterLonModHalfPi,
        -sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
    double latRad = asin(zDxDyEq1); 
    double deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi + deltaLon <= PI_OVER_FOUR) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad + deltaLon, latRad);
    } else { // Look at the neighbour base cell, comute north east
      zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          -sinConeCenterLat, true, coneCenterLonModHalfPi - HALF_PI,
          -sinConeCenterLat, twoSineOfHalfConeRadius, false, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad + deltaLon, latRad) : - 1L;
    }
    // North-west
    zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
        zStart, false, coneCenterLonModHalfPi,
        -sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
    latRad = asin(zDxDyEq1); 
    deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
        coneCenterLatRad - latRad, zDxDyEq1);
    if (isFinite(deltaLon) && coneCenterLonModHalfPi - deltaLon >= 0) {
      result[resultSize++] = hashComputer.hash(coneCenterLonRad - deltaLon, latRad);
    } else { // Look at the neighbour base cell, comute south west
      zDxDyEq1 = -NewtonMethod.newtonSolveNorthPolarCapZone(
          -sinConeCenterLat, false, HALF_PI + coneCenterLonModHalfPi,
          -sinConeCenterLat, twoSineOfHalfConeRadius, true, relativePrecision, nIterMax);
      latRad = asin(zDxDyEq1);
      deltaLon = angDistComputer.coneDeltaLon(squareOfsinOfHalfR, cosConeCenterLat,
          coneCenterLatRad - latRad, zDxDyEq1);
      result[resultSize++] = isFinite(deltaLon) ? hashComputer.hash(coneCenterLonRad - deltaLon, latRad): -1L;
    }
    return resultSize;
  }

}
