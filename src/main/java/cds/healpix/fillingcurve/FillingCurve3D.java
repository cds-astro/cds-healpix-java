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

package cds.healpix.fillingcurve;

/**
 * Interface defining a 3-dimensional space-filling curve.
 * It provides method to transform 3d-coordinates in a single index in the curve.
 * The index value is here called hash value since the intput parameters are discritized
 * coordinates, and information is lost during the discretization.  
 * 
 * @author F.-X. Pineau
 *
 */
public interface FillingCurve3D {

  /**
   * Transforms coordinates in a discretized 2d-space into a space
   * filling curve index. This can be thought of as computing a hash
   * value from 3d-coordinates.
   * @param i discretized coordinate along the horizontal axis
   * @param j discretized coordinate along the vertical axis
   * @param k discretized coordinate along the depth axis
   * @return the space filling curve index (or hash value) associated
   *     to the given discretized 3d-coordinates.
   */
  long ijk2hash(int i, int j, int k);

  /**
   * Special case of {@link #ijk2hash(int, int, int)} in which the discretized coordinate along
   * the vertical axis equals zero.
   * @param i discretized coordinate along the horizontal axis
   * @return the space filling curve index (or hash value) associated to the given discretized
   *     horizontal coordinate, assuming the discretized vertical coordinate equal zero.
   */
  long i002hash(int i);

  /**
   * Transforms the given space filling curve index (or hash value) into
   * a single value from which it is straightforward to extract the
   * associated 3d-coordinates using methods {@link #ijk2i(long)}, {@link #ijk2j(long)} and
   * {@link #ijk2k(long)}.
   * @param hash the space filling curve index (or hash value)
   * @return a single value from which it is straightforward to extract the associated
   *     3d-coordinates using methods {@link #ijk2i(long)}, {@link #ijk2j(long)} and
   * {@link #ijk2k(long)}.
   */
  long hash2ijk(long hash);

  /**
   * Special case of {@link #hash2ijk(long)} in which the discretized coordinate along
   * the vertical axis equals zero.
   * @param hash the space filling curve index (or hash value)
   * @return a single value, knowing the vertical coordinate equals zero, from which it is
   *     straightforward to extract the associated horizontal coordinate using method
   *     {@link #ijk2i(long)}.
   */
  long hash2i00(long hash);

  /**
   * Extract the discretized horizontal coordinates from the result of the method
   * {@link #hash2ijk(long)}.
   * @param ijk result of the method {@link #hash2ijk(long)}
   * @return the discretized horizontal coordinate stored in the given parameter {@code ijk}
   */
  int ijk2i(long ijk);

  /**
   * Extract the discretized vertical coordinates from the result of the method
   * {@link #hash2ijk(long)}.
   * @param ijk result of the method {@link #hash2ijk(long)}
   * @return the discretized horizontal coordinate stored in the given parameter {@code ijk}
   */
  int ijk2j(long ijk);

  /**
   * Extract the discretized depth coordinates from the result of the method
   * {@link #hash2ijk(long)}.
   * @param ijk result of the method {@link #hash2ijk(long)}
   * @return the discretized horizontal coordinate stored in the given parameter {@code ijk}
   */
  int ijk2k(long ijk);

}
