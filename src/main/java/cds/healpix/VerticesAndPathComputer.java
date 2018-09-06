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

/**
 * Interface defining methods used to compute the position (on the unit sphere) of HEALPix cells
 * center and vertices, and paths along HEALPix cells sides and edge. 
 * 
 * @author F.-X. Pineau
 *
 */
public interface VerticesAndPathComputer extends HierarchyItem {

  /** {@link EnumSet} containing all elements of the {@link Cardinal} enum.
   * Defined for users comodity. */
  public static final EnumSet<Cardinal> ALL_CARDINAL_POINTS = EnumSet.allOf(Cardinal.class);
  
  /** Index of the longitude in the array containing coordinate on the unit sphere. */
  static final int LON_INDEX = Projection.LON_INDEX;
  /** Index of the lattitude in the array containing coordinate on the unit sphere. */
  static final int LAT_INDEX = Projection.LAT_INDEX;
  
  /**
   * Compute the position on the unit sphere of the center (in the Euclidean projection plane)
   * of the cell associated to the given hash value.
   * @param hash hash value of the cell we look for the unprojected center
   * @return the unprojected position (on the unit sphere) of the center of the cell in the
   *         Euclidean plane. The lon and lat coordinate are stored in the returned array at indices
   *         {@link #LON_INDEX} and {@link #LAT_INDEX} respectively.
   *         Lat in [-pi/2, pi/2] radians, lon is in [0, 2pi] radians.
   *          
   */
  double[] center(long hash);

  /**
   * See {@link #center(long)}, except that the result is stored in the given array.
   * @param hash hash value of the cell we look for the unprojected center.
   * @param resultLonLat array used to store the result. Must be of size &gt;= 2.
   */
  void center(long hash, double[] resultLonLat);
  

  /**
   * Compute the position of an HEALPix cell vertex on the unit sphere.
   * @param hash hash value of the cell we look for the given vertex.
   * @param cardinalPoint location of the vertex with respect to the cell center
   * @return the position (on the unit sphere) of the vertex located at the given cardinal direction
   *         from the center of the given cell.
   *         The lon and lat coordinate are stored in the returned array at indices
   *         {@link #LON_INDEX} and {@link #LAT_INDEX} respectively.
   *         Lat is in [-pi/2, pi/2] radians, lon is in [0, 2pi] radians.
   */
  double[] vertex(long hash, Cardinal cardinalPoint);

  /**
   * See {@link #vertex(long, Cardinal)}, except that the result is stored in the given array.
   * @param hash hash value of the cell we look for the given vertex.
   * @param cardinalPoint location of the vertex with respect to the cell center
   * @param resultLonLat array used to store the result. Must be of size &gt;= 2.
   */
  void vertex(long hash, Cardinal cardinalPoint, double[] resultLonLat);
  
  /**
   * Returns the vertices located at the given cardinal points.
   * If you want the full 4 vertices, simply use {@code EnumSet.allOf(Cardinal)}.
   * For West and East vertices, use {@code EnumSet.of(Cardinal.W, Cardinal.E)}. 
   * @param hash hash value of the cell we look for the given vertices.
   * @param cardinalPoints locations of the vertices we look for the positions
   * @return the positions (on the unit sphere) of the vertices located at the given cardinal
   *         directions from the center of the given cell.
   *         The lon and lat coordinates of each vertex are stored in the arrays at indices
   *         {@link #LON_INDEX} and {@link #LAT_INDEX} respectively.
   *         Lat is in [-pi/2, pi/2] radians, lon is in [0, 2pi] radians.
   */
  EnumMap<Cardinal, double[]> vertices(long hash, EnumSet<Cardinal> cardinalPoints);
  
  /**
   * See {@link #vertices(long, EnumSet)}. The difference is that the user provides a
   * pre-set Map. The structure of the Map is not modified, but the coordinates of the values
   * (array values) are overwritten. An error will be thrown if a value is null or contains
   * less than two elements.
   * @param hash hash value of the cell we look for the given vertices.
   * @param cardinalPoints the map to be modified y
   */
  void vertices(long hash, final EnumMap<Cardinal, double[]> cardinalPoints);
  
  /**
   * Compute points on a given side of a given HEALPix cell on the unit sphere.
   * @param hash hash value of the cell we look for side path on the unit sphere.
   * @param fromVertex direction (from the cell center) of the path starting vertex
   * @param toVertex   direction (from the cell center) of the path ending vertex
   * @param isToVertexIncluded if {@code false}, the result contains {@code nSegments} points and do
   *                           no include the result ending vertex. Else the result contains
   *                           {@code nSegments + 1} points.
   * @param nSegments number of segments in the path from the starting vertex to the ending vertex
   * @return a list of points on the given side of the given HEALPix cell on the unit sphere.
   */
  double[][] pathAlongCellSide(long hash, Cardinal fromVertex, Cardinal toVertex,
      boolean isToVertexIncluded, int nSegments);
  
  /**
   * See {@link #pathAlongCellSide(long, Cardinal, Cardinal, boolean, int)}. The difference is that
   * the user provides a list of points whose coordinates are going to be overwritten.
   * An error will be thrown if the array (of array) is not large enough (i.e. its is smaller than 
   * nSegments or nSegments + 1) or if one of the array is null or contains less than two elements.
   * @param hash hash value of the cell we look for side path on the unit sphere.
   * @param fromVertex direction (from the cell center) of the path starting vertex
   * @param toVertex   direction (from the cell center) of the path ending vertex
   * @param isToVertexIncluded if {@code false}, the result contains {@code nSegments} points and do
   *                           no include the result ending vertex. Else the result contains
   *                           {@code nSegments + 1} points.
   * @param nSegments number of segments in the path from the starting vertex to the ending vertex
   * @param pathPoints object used to store the result.
   */
  void pathAlongCellSide(long hash, Cardinal fromVertex, Cardinal toVertex, 
      boolean isToVertexIncluded, int nSegments, double[][] pathPoints);
  

  /**
   * 
   * @param hash hash value of the cell we look for the edge path on the unit sphere.
   * @param startingVertex direction (from the cell center) of the path starting vertex
   * @param clockwiseDirection if {@code true}, result points are stored in a clockwise direction
   *                           order, else they are stored in counter-clockwise direction order.
   * @param nSegmentsBySide number of segments in each each side is divided. Hence, the total number
   *                        of points in the result is {@code 4 x nSegmentsBySide}.
   * @return a list of points on the given edge of the given HEALPix cell on the unit sphere,
   *        clockwise or counter-clockwise ordered.
   */
  double[][] pathAlongCellEdge(long hash, Cardinal startingVertex, boolean clockwiseDirection,
      int nSegmentsBySide);
  
  /**
   * See {@link #pathAlongCellEdge(long, Cardinal, boolean, int)}. The difference is that the user
   * provides a list of points whose coordinates are going to be overwritten.
   * An error will be thrown if the array (of array) is not large enough (i.e. its is smaller than 
   * 4 *nSegments) or if one of the array is null or contains less than two elements.
   * @param hash hash hash value of the cell we look for the edge path on the unit sphere.
   * @param startingVertex direction (from the cell center) of the path starting vertex
   * @param clockwiseDirection if {@code true}, result points are stored in a clockwise direction
   *                           order, else they are stored in counter-clockwise direction order.
   * @param nSegmentsBySide number of segments in each each side is divided. Hence, the total number
   *                        of points in the result is {@code 4 x nSegmentsBySide}.
   * @param pathPoints object used to store the result.
   */
  void pathAlongCellEdge(long hash, Cardinal startingVertex, boolean clockwiseDirection,
      int nSegmentsBySide, double[][] pathPoints);
  
}
