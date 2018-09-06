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

final class BaseHashes {

  private static final BaseHash[] BASE_HASHES = new BaseHash[12];
  static {
    BASE_HASHES[ 0] = new BaseHashNorthPolarCap(0);
    BASE_HASHES[ 1] = new BaseHashNorthPolarCap(1);
    BASE_HASHES[ 2] = new BaseHashNorthPolarCap(2);
    BASE_HASHES[ 3] = new BaseHashNorthPolarCap(3);
    BASE_HASHES[ 4] = new BaseHashEquatorial(4);
    BASE_HASHES[ 5] = new BaseHashEquatorial(5);
    BASE_HASHES[ 6] = new BaseHashEquatorial(6);
    BASE_HASHES[ 7] = new BaseHashEquatorial(7);
    BASE_HASHES[ 8] = new BaseHashSouthPolarCap(8);
    BASE_HASHES[ 9] = new BaseHashSouthPolarCap(9);
    BASE_HASHES[10] = new BaseHashSouthPolarCap(10);
    BASE_HASHES[11] = new BaseHashSouthPolarCap(11);
  }
  
  /** Prevents from outside instantiation. */
  private BaseHashes() { }
  
  static final BaseHash get(int d0h) {
    return BASE_HASHES[d0h];
  }

}
