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

import static cds.healpix.CompassPoint.MainWind.*;

import java.util.EnumSet;

import cds.healpix.CompassPoint.MainWind;

final class BaseHashEquatorial extends BaseHash {

  private static final EnumSet<MainWind> noNeighbour = EnumSet.of(N, S);
  
  protected BaseHashEquatorial(int baseHash) {
    super(baseHash);
  }

  @Override
  protected void fillNeighbours(int d0h, int[] neighbours) {
    assert 4 <= d0h && d0h < 8;
    final int hModulo4 = d0h & 3;
    int hp1Modulo4 = (hModulo4 + 1) & 3;
    neighbours[SE.getIndex()] = hModulo4 + 8;
    neighbours[ E.getIndex()] = hp1Modulo4 + 4;
    neighbours[NE.getIndex()] = hModulo4;
    neighbours[ S.getIndex()] = -1; // No neighbour in South
    neighbours[ C.getIndex()] = d0h;
    neighbours[ N.getIndex()] = -1; // No neighbour in North
    int hm1Modulo4 = (hModulo4 - 1) & 3;
    neighbours[NW.getIndex()] = hm1Modulo4;
    neighbours[ W.getIndex()] = hm1Modulo4 + 4;
    neighbours[SW.getIndex()] = hm1Modulo4 + 8;
  }

  @Override
  MainWind getDirectionFromNeighbour(final MainWind neighbourMainWind) {
    if (noNeighbour.contains(neighbourMainWind)) {
      throw new IllegalArgumentException("No neighbour in direction " + neighbourMainWind);
    }
    return neighbourMainWind.getOppositeDirection();
  }
  
  @Override
  int pickRightIndexOnNeighbourSouthToEastAxis(MainWind neighbourDirection,
      int iAxisSE, int iAxisSW, int nsideMinus1) {
    return neighbourDirection.pickRightIntValue(
            iAxisSE,      -1,       0,
        nsideMinus1, iAxisSE,       0,
        nsideMinus1,      -1, iAxisSE);
  }

  @Override
  long pickRightIndexBitsOnNeighbourSouthToEastAxis(MainWind neighbourDirection,
      long iAxisSEBits, long jAxisSWBits, long nsideMinus1Bits) {
    return neighbourDirection.pickRightLongValue(
            iAxisSEBits,          -1,           0,
        nsideMinus1Bits, iAxisSEBits,           0,
        nsideMinus1Bits,          -1, iAxisSEBits);
  }

  @Override
  int pickRightIndexOnNeighbourSouthToWestAxis(MainWind neighbourDirection,
      int iAxisSE, int iAxisSW, int nsideMinus1) {
    return neighbourDirection.pickRightIntValue(
              0,      -1,     iAxisSW,
              0, iAxisSW, nsideMinus1,
        iAxisSW,      -1, nsideMinus1);
  }

  @Override
  long pickRightIndexSWBits(MainWind neighbourDirection,
      long iAxisSEBits, long jAxisSWBits, long nsideMinus1Bits) {
    return neighbourDirection.pickRightLongValue(
                   0,          -1,    jAxisSWBits,
                   0, jAxisSWBits, nsideMinus1Bits,
         jAxisSWBits,          -1, nsideMinus1Bits);
  }

}
