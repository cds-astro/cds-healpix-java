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

import cds.healpix.CompassPoint.Cardinal;

interface ConeOrdinalHashComputer {

  /**
   * 
   * @param coneCenterLonRad
   * @param coneCenterLatRad
   * @param coneRadiusRad
   * @param hashComputer
   * @param angDistComputer
   * @param relativePrecision
   * @param nIterMax
   * @param hashCenterAtSmallestDepth
   * @param verticesComputerAtSmallestDepth
   * @param vertices passsed to avoid creating a new objet at each call. WARNING: its content is overwritten!
   * @param cosConeCenterLat
   * @param sinConeCenterLat
   * @param twoSineOfHalfConeRadius
   * @param squareOfsinOfHalfR
   * @param result
   * @return
   */
  int computeOrdinalHash(
      // Minimum parameters required
      double coneCenterLonRad, double coneCenterLatRad,
      double coneRadiusRad, HashComputer hashComputer, AngularDistanceComputer angDistComputer,
      // Alog params
      double relativePrecision, int nIterMax,
      // Pre-computed quantities
      final long hashCenterAtSmallestDepth,
      final VerticesAndPathComputer verticesComputerAtSmallestDepth, 
      final EnumMap<Cardinal, double[]> vertices,
      final double cosConeCenterLat, final double sinConeCenterLat,
      final double twoSineOfHalfConeRadius, final double squareOfsinOfHalfR,
      // Store the result
      long[] result);
  
}
