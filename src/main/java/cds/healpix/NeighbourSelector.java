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

import java.util.EnumSet;

import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.CompassPoint.MainWind;
import cds.healpix.CompassPoint.Ordinal;

/**
 * 
 * In this interface methods, we do not used an {@code EnumMap} of {@link CompassPoint.MainWind} but a
 * {@link NeighbourList} because:
 * 1 - we want to avoid autoboxing (a map storing Long and not long).
 * 2 - we want to transparently deal with cases in which there is no neighbour in a given direction. 
 * 
 * @author F.-X. Pineau
 *
 */
public interface NeighbourSelector {

  static final EnumSet<MainWind> ALL_MAIN_WINDS = EnumSet.allOf(MainWind.class);
  
  /**
   * The hash value of the neighbour of the cell of given hash, in the given deirection.
   * @param hash ash value of the cell we are looking for the neighbour.
   * @param direction direction of the neighbour we are looking for.
   * @return -1 if there is no neighbour in the given direction (hashes located at corners of
   * base hases).
   */
  long neighbour(long hash, MainWind direction);


  /**
   * Returns the list of the hash of the cells surrounding the cell defined
   * by the given hash. The number of surrounding cells can be 8 (for all
   * cells except the cells located at the corners of the 12 depth 0 cells), 
   * 7 (for the west and east corners of the polar caps, i.e. depth 0 cells
   * number 0, 1, 2, 4, 8, 9, 10 and 11, or for the north and south corners
   * of the equatorial regions, i.e. depth 0 cells number 4, 5, 6 and 7), or
   * 6 for depth 0 pixels. 
   * @param hash hash value of the cell we are looking for the neighbours.
   * @return the list of neighbours hashes
   */
  NeighbourList neighbours(long hash);

  /**
   * Equivalent of {@link #neighbours(long)} but passing in argument the 
   * list object to be filled. 
   * WARNING: the content of the provided list is overwritten, and be sure
   * the list is large enough (at least of size 8).
   * And the value -1 is returned if the given hash has no neighbour at the given main wind direction.
   * @param hash the hash code we want the neighbours
   * @param result which contains the list of neighbours from index 0 to the number of neighbours
   */
  void neighbours(long hash, NeighbourList result);

  /**
   * Idem as {@link #neighbours(long, NeighbourList)} except that the result is put in a 
   * simple {@link FlatHashList}.
   * WARNING: the content of the provided list is overwritten, and be sure
   * the list is large enough (at least of size 8).
   * @param hash the hash code we want the neighbours
   * @param result  which contains the list of neighbours from index 0 to the number of neighbours
   */
  void neighbours(long hash, FlatHashList result);
  
  /**
   * Equivalent of {@link #neighbours(long, FlatHashList)} but with the list of direction of the
   * wanted neighbours.
   * @param hash hash value of the cell we are looking for the neighbours.
   * @param directions the directions of the neighbours we are looking for.
   * @return the list of the hash of the neighbours of the cell having the given hash.
   */
  NeighbourList neighbours(long hash, EnumSet<MainWind> directions);

  /**
   * Equivalent of {@link #neighbours(long, FlatHashList)} but with the list of direction of the
   * wanted neighbours.
   * And the value -1 is returned if the given hash has no neighbour at the given main wind direction.
   * @param hash hash value of the cell we are looking for the neighbours.
   * @param directions the directions of the neighbours we are looking for.
   * @param result the list holding the result.
   */
  void neighbours(long hash, EnumSet<MainWind> directions, NeighbourList result);

  
  
  FlatHashList internalEdges(long hash, int deltaDepth);
  void internalEdges(long hash, int deltaDepth, FlatHashList result);
          
  
  FlatHashList sortedInternalEdges(long hash, int deltaDepth);
  void sortedInternalEdges(long hash, int deltaDepth, FlatHashList result);
          
          
  FlatHashList sortedInternalEdge(long hash, int deltaDepth, Ordinal direction);
  void sortedInternalEdge(long hash, int deltaDepth, Ordinal direction, FlatHashList result);
  
  FlatHashList sortedInternalEdgeNE(long hash, int deltaDepth);
  void sortedInternalEdgeNE(long hash, int deltaDepth, FlatHashList result);
  
  FlatHashList sortedInternalEdgeNW(long hash, int deltaDepth);
  void sortedInternalEdgeNW(long hash, int deltaDepth, FlatHashList result);
  
  FlatHashList sortedInternalEdgeSE(long hash, int deltaDepth);
  void sortedInternalEdgeSE(long hash, int deltaDepth, FlatHashList result);
  
  FlatHashList sortedInternalEdgeSW(long hash, int deltaDepth);
  void sortedInternalEdgeSW(long hash, int deltaDepth, FlatHashList result);
  
  long internalCorner(long hash, int deltaDepth, Cardinal direction);
  long internalCornerN(long hash, int deltaDepth);
  long internalCornerS(long hash, int deltaDepth);
  long internalCornerE(long hash, int deltaDepth);
  long internalCornerW(long hash, int deltaDepth);

  FlatHashList externalEdges(long hash, int deltaDepth);
  void externalEdges(long hash, int deltaDepth,  final FlatHashList result);
  
  FlatHashList sortedExternalEdges(long hash, int deltaDepth);
  void sortedExternalEdges(long hash, int deltaDepth,  final FlatHashList result);

}
