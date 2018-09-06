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
 * Compute an HEALPix hash value additionaly taking a 3rd dimension (the auxiliary axis).
 * 
 * @author F.-X. Pineau
 *
 */
public interface HashComputerWithAux extends HierarchyItem {
  /**
   * Returns the depth of the hash computed by {@link #hash(double, double, double)} method.
   * @return the depth of the hash computed by {@link #hash(double, double, double)} method.
   */
  @Override
  int depth();

  /**
   * Returns the HEALPix hash value of the given coordinate at this object depth.
   * WARNING: depending on the implementation, this method may or may not be thread-safe.
   * @param lonRad longitude in radians, must support reasonably large positive and negative values
   *         producing accurate results with a naive range reduction like modulo 2*pi
   *         (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
   * @param latRad latitude in [-pi/2, pi/2] radians
   * @param auxValue value on the 3rd (the auxiliary) axis.
   * @return the hash value associated to the given coordinate, at this object depth.
   */
  long hash(final double lonRad, final double latRad, final double auxValue);
}
