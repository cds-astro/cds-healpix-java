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

import static cds.healpix.CompassPoint.MainWind.C;
import static cds.healpix.CompassPoint.MainWind.E;
import static cds.healpix.CompassPoint.MainWind.N;
import static cds.healpix.CompassPoint.MainWind.NE;
import static cds.healpix.CompassPoint.MainWind.NW;
import static cds.healpix.CompassPoint.MainWind.S;
import static cds.healpix.CompassPoint.MainWind.SE;
import static cds.healpix.CompassPoint.MainWind.SW;
import static cds.healpix.CompassPoint.MainWind.W;

import java.util.EnumMap;
import java.util.EnumSet;

import cds.healpix.CompassPoint.MainWind;

final class BaseHashSouthPolarCap extends BaseHash {

  private static final EnumSet<MainWind> noNeighbour = EnumSet.of(E, W);
  private static final EnumMap<MainWind, MainWind> dirFromNeighbour = new EnumMap<MainWind, MainWind>(MainWind.class);
  static {
    dirFromNeighbour.put(NE, NE.getOppositeDirection());
    dirFromNeighbour.put(N,   N.getOppositeDirection());
    dirFromNeighbour.put(NW, NW.getOppositeDirection());
    dirFromNeighbour.put(C, C.getOppositeDirection());    
    dirFromNeighbour.put(SE, SW);
    dirFromNeighbour.put(S,   S);
    dirFromNeighbour.put(SW, SE);
  }
  
  protected BaseHashSouthPolarCap(int baseHash) {
    super(baseHash);
  }

  @Override
  MainWind getDirectionFromNeighbour(final MainWind neighbourMainWind) {
    if (noNeighbour.contains(neighbourMainWind)) {
      throw new IllegalArgumentException("No neighbour in direction " + neighbourMainWind);
    }
    return dirFromNeighbour.get(neighbourMainWind);
  }
  
  @Override
  protected void fillNeighbours(int d0h, int[] neighbours) {
    assert 8 <= d0h && d0h < 12;
    final int hModulo4 = d0h & 3;
    int hp1Modulo4 = (hModulo4 + 1) & 3;
    neighbours[SE.getIndex()] = hp1Modulo4 + 8;
    neighbours[ E.getIndex()] = -1; // No neighbour  in East
    neighbours[NE.getIndex()] = hp1Modulo4 + 4;
    int hp2Modulo4 = (hModulo4 + 2) & 3;
    neighbours[ S.getIndex()] = 8 + hp2Modulo4;
    neighbours[ C.getIndex()] = d0h;
    neighbours[ N.getIndex()] = hModulo4;
    int hm1Modulo4 = (hModulo4 - 1) & 3;
    neighbours[NW.getIndex()] = hModulo4 + 4;
    neighbours[ W.getIndex()] = -1; // No neighbour  in West
    neighbours[SW.getIndex()] = hm1Modulo4 + 8;
  }

  @Override
  int pickRightIndexOnNeighbourSouthToEastAxis(MainWind neighbourDirection, 
      int iAxisSE, int jAxisSW, int nsideMinus1) {
    return neighbourDirection.pickRightIntValue(
        iAxisSE,       0,  0,
             -1, iAxisSE, -1,
        jAxisSW,       0,  0);
  }

  @Override
  long pickRightIndexBitsOnNeighbourSouthToEastAxis(MainWind neighbourDirection,
      long iAxisSEBits, long jAxisSWBits, long nsideMinus1Bits) {
    return neighbourDirection.pickRightLongValue(
      iAxisSEBits,           0,  0,
               -1, iAxisSEBits, -1,
      jAxisSWBits,           0,  0);
  }

  @Override
  int pickRightIndexOnNeighbourSouthToWestAxis(MainWind neighbourDirection,
      int iAxisSE, int jAxisSW, int nsideMinus1) {
    return neighbourDirection.pickRightIntValue(
         0,       0, jAxisSW,
        -1, jAxisSW,      -1,
         0,       0, iAxisSE);
  }

  @Override
  long pickRightIndexSWBits(MainWind neighbourDirection,
      long iAxisSEBits, long jAxisSWBits, long nsideMinus1Bits) {
    return neighbourDirection.pickRightLongValue(
         0,           0, jAxisSWBits,
        -1, jAxisSWBits,          -1,
         0,           0, iAxisSEBits);
  }
  
}
