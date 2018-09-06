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

/**
 * The idea of this interface is to avoid making multiple time the same operations (like selecting
 * the optimal starting depth) in case of fixed radius cross-match.
 * 
 * @author F.-X. Pineau
 *
 */
public interface HealpixNestedFixedRadiusConeComputer {

  /**
   * Returns the radius of the cones, in radians.
   * @return he radius of the cones, in radians.
   */
  double getRadius();
  
  /**
   * MOC of the cells having a part of their surface area in common with the given cone
   * @param coneCenterLonRad longitude of the center of the cone, in radians
   * @param coneCenterLatRad latitude of the center of the cone, in radians
   * @return the resulting MOC.
   */
  HealpixNestedBMOC overlappingCells(double coneCenterLonRad, double coneCenterLatRad);
  
  /**
   * MOC of the cells which centers lie inside the given cone.
   * @param coneCenterLonRad longitude of the center of the cone, in radians
   * @param coneCenterLatRad latitude of the center of the cone, in radians
   * @return the resulting MOC.
   */
  HealpixNestedBMOC overlappingCenters(double coneCenterLonRad, double coneCenterLatRad);

  
  /**
   * To obtain new instances in case we want to use multi-threading,
   * since an object is possibly not thread-safe.
   * If the objectis thread-safe, the method can simply return {@code this}.
   * @return a new instance of {@link HealpixNestedFixedRadiusConeComputer}.
   */
  HealpixNestedFixedRadiusConeComputer newComputer();
  
}
