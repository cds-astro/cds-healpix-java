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

final class BaseHashNorthPolarCap extends BaseHash {

  private static final EnumSet<MainWind> noNeighbour = EnumSet.of(E, W);
  private static final EnumMap<MainWind, MainWind> dirFromNeighbour = new EnumMap<MainWind, MainWind>(MainWind.class);
  static {
    dirFromNeighbour.put(NE, NW);
    dirFromNeighbour.put(N,   N);
    dirFromNeighbour.put(NW, NE);
    dirFromNeighbour.put(C, C.getOppositeDirection());
    dirFromNeighbour.put(SE, SE.getOppositeDirection());
    dirFromNeighbour.put(S,   S.getOppositeDirection());
    dirFromNeighbour.put(SW, SW.getOppositeDirection());
  }
  
  protected BaseHashNorthPolarCap(int baseHash) {
    super(baseHash);
  }

  @Override
  protected void fillNeighbours(int d0h, int[] neighbours) {
    assert 0 <= d0h && d0h < 4;
    final int hModulo4 = d0h & 3; 
    int hp1Modulo4 = (d0h + 1) & 3;
    neighbours[SE.getIndex()] = hp1Modulo4 + 4;
    neighbours[ E.getIndex()] = -1; // No eastern neighbour
    neighbours[NE.getIndex()] = hp1Modulo4;
    int hp2Modulo4 = (d0h + 2) & 3;
    neighbours[ S.getIndex()] = d0h + 8;
    neighbours[ C.getIndex()] = d0h;
    neighbours[ N.getIndex()] = hp2Modulo4;
    int hm1Modulo4 = (d0h - 1) & 3;
    neighbours[NW.getIndex()] = hm1Modulo4;
    neighbours[ W.getIndex()] = -1; // No western neighbour
    neighbours[SW.getIndex()] = hModulo4 + 4;
  }

  @Override
  MainWind getDirectionFromNeighbour(final MainWind neighbourMainWind) {
    if (noNeighbour.contains(neighbourMainWind)) {
      throw new IllegalArgumentException("No neighbour in direction " + neighbourMainWind);
    }
    return dirFromNeighbour.get(neighbourMainWind);
  }
  
  @Override
  int pickRightIndexOnNeighbourSouthToEastAxis(MainWind neighbourDirection,
      int iAxisSE, int iAxisSW, int nsideMinus1) {
    return neighbourDirection.pickRightIntValue(
        nsideMinus1, nsideMinus1,  iAxisSW,
                 -1,     iAxisSE,       -1,
        nsideMinus1,  nsideMinus1, iAxisSE);
  }

  @Override
  long pickRightIndexBitsOnNeighbourSouthToEastAxis(MainWind neighbourDirection,
      long iAxisSEBits, long jAxiSWBits, long nsideMinus1Bits) {
    return neighbourDirection.pickRightLongValue(
        nsideMinus1Bits, nsideMinus1Bits, jAxiSWBits,
                     -1,     iAxisSEBits,          -1,
        nsideMinus1Bits, nsideMinus1Bits, iAxisSEBits);
  }

  @Override
  int pickRightIndexOnNeighbourSouthToWestAxis(MainWind neighbourDirection,
      int iAxisSE, int iAxisSW, int nsideMinus1) {
    return neighbourDirection.pickRightIntValue(
        iAxisSE, nsideMinus1, nsideMinus1,
            -1,      iAxisSW,          -1,
        iAxisSW, nsideMinus1, nsideMinus1);
  }

  @Override
  long pickRightIndexSWBits(MainWind neighbourDirection,
      long iAxisSEBits, long jAxisSWBits, long nsideMinus1Bits) {
    return neighbourDirection.pickRightLongValue(
        iAxisSEBits, nsideMinus1Bits, nsideMinus1Bits,
           -1,           jAxisSWBits,              -1,
        jAxisSWBits, nsideMinus1Bits, nsideMinus1Bits);
  }

}
