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
 * Define a projection of spherical coordinates (in the unit sphere) to the Euclidean plane, together
 * with the reverse operation.
 * Remark: we use arrays instead of objects because HEALPix is a low level library focused on 
 *         performances. Accessing an element in an array is faster than accessing a class attribute. 
 * 
 * @author F.-X. Pineau
 *
 */
public interface Projection {

  /** Index of the longitude in the array containing the result of an unproject method. */
  static final int LON_INDEX = 0;
  /** Index of the lattitude in the array containing the result of an unproject method. */
  static final int LAT_INDEX = 1;
  
  /** Index of the x coordinate in the array containing the result of a project method. */
  static final int X_INDEX = 0;
  /** Index of the y coordinate in the array containing the result of a project method. */
  static final int Y_INDEX = 1;
  
  /**
   * Project the given spherical coordinates into the Euclidean plane.
   * @param lonRad longitude in radians (the accepted value range is implementation dependent)
   * @param latRad latitude in [-pi/2, pi/2] radians
   * @return the projection of the given spherical coordinate into the Euclidean plane.
   *         The x and y coordinate are stored in the returned array at indices
   *         {@link #X_INDEX} and {@link #Y_INDEX} respectively.
   *         The range of possible value for x and y is implementation dependent.
   */
  double[] project(double lonRad, double latRad);
  
  /**
   * See {@link #project(double, double)} with the result stored in the given array.
   * @param lonRad see {@link #project(double, double)}
   * @param latRad see {@link #project(double, double)}
   * @param resultXY array used to store the result. Must be of size &gt;= 2.
   */
  void project(double lonRad, double latRad, double[] resultXY);
  
  /**
   * Reverse projection: we look for spherical coordinates from their projected coordinates
   * ({@link #project(double, double)}.  
   * @param x the x coordinate of the projected spherical point we are looking for,
   *          the accepted value range is implementation dependent
   * @param y the y coordinate of the projected spherical point we are looking for,
   *          must be in [-2, 2].
   * @return the spherical coordinates leading to the given projection coordinates.
   *         The lon and lat coordinate are stored in the returned array at indices
   *         {@link #LON_INDEX} and {@link #LAT_INDEX} respectively.
   *         Lat is in [-pi/2, pi/2] radians, lon is also in radians but it possible value range
   *         is implementation dependent).
   */
  double[] unproject(double x, double y);
  
  /**
   * See {@link #unproject(double, double)} with the result stored in the given array.
   * @param x {@link #unproject(double, double)}
   * @param y {@link #unproject(double, double)}
   * @param resultLonLat array used to store the result. Must be of size &gt;= 2.
   */
  void unproject(double x, double y, double[] resultLonLat);

}
