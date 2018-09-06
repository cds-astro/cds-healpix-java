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

import cds.healpix.CompassPoint.MainWind;

/**
 * Represent one of the 12 hash of depth zero (often called d0h in the code).
 * 
 * @author F.-X. Pineau
 *
 */
abstract class BaseHash {

  private final int d0h;

  private final int[] neighbours;

  protected BaseHash(final int baseHash) {
    this.d0h = baseHash;
    this.neighbours = new int[MainWind.size()];
    this.fillNeighbours(this.d0h, this.neighbours);
  }

  protected abstract void fillNeighbours(int d0h, int[] neighbours);

  final int getValue() {
    return this.d0h;
  }

  // return -1 if there is no neighbour at the given main wind direction.
  final int getNeighbour(MainWind mainWind) {
    return this.neighbours[mainWind.getIndex()];
  }

  /**
   * 
   * @param neighbourMainWind
   * @return the direction of this base cell from its neighbour located at the given  direction.
   */
  abstract MainWind getDirectionFromNeighbour(MainWind neighbourMainWind);
  
  /**
   * Method used when looking for neighbours of hashes which are on the border (including corners)
   * of a base hash. 
   * @param neighbourDirection direction of the neighbour base cell
   * @param indexOnAxisS2E index in the current base hash, must be in [-1, nside]
   * @param indexOnAxisS2W index in the current base hash, must be in [-1, nside]
   * @param nsideMinus1 nside -1, i.e the max value of an index on the South-East or South-West axis
   * @return the index on the neighbour South-East axis of the hash of given position in the current
   * base hash
   */
  abstract int pickRightIndexOnNeighbourSouthToEastAxis(MainWind neighbourDirection,
      int indexOnAxisS2E, int indexOnAxisS2W, int nsideMinus1);
  
  /**
   * Method used when looking for neighbours of hashes which are on the border (including corners)
   * of a base hash. 
   * @param neighbourDirection direction of the neighbour base cell
   * @param inedxOnAxisS2E index in the current base hash, must be in [-1, nside]
   * @param indexOnAxisS2W index in the current base hash, must be in [-1, nside]
   * @param nsideMinus1 nside -1, i.e the max value of an index on the South-East or South-West axis
   * @return the index on the neighbour South-West axis of the hash of given position in the current
   * base hash
   */
  abstract int pickRightIndexOnNeighbourSouthToWestAxis(MainWind neighbourDirection,
      int inedxOnAxisS2E, int indexOnAxisS2W, int nsideMinus1);
  
  abstract long pickRightIndexBitsOnNeighbourSouthToEastAxis(MainWind neighbourDirection, 
      long indexBitsOnAxisS2E, long indexBitsOnAxisS2W, long nsideMinus1Bits);
  
  abstract long pickRightIndexSWBits(MainWind neighbourDirection,
      long indexBitsOnAxisS2E, long indexBitsOnAxisS2W, long nsideMinus1Bits);

}
