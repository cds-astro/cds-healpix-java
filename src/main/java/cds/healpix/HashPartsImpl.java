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
 * Data-structure storing the base hash value (depth 0 hash value) of the hash it represents, 
 * together with the (i, j) coordinates of the hash in its base hash sub-divided "depth" time.
 * The "depth" is not stored explicitly in the data-structure.
 * 
 * @author F.-X. Pineau
 *
 */
final class HashPartsImpl implements HashParts {
  
  final int d0h;
  final int iInD0h;
  final int jInD0h;
  
  HashPartsImpl(final int baseHash, final int iInBasePixel, final int jInBasePixel) {
    this.d0h = baseHash;
    this.iInD0h = iInBasePixel;
    this.jInD0h = jInBasePixel;
  }

  @Override
  public int baseCellHash() {
    return this.d0h;
  }

  @Override
  public int iInBaseCell() {
    return this.iInD0h;
  }

  @Override
  public int jInBaseCell() {
    return this.jInD0h;
  }
  
  @Override
  public String toString() {
    return "[d0h: " + this.d0h + ", i:" + this.iInD0h + ", j: " + this.jInD0h + "]";
  }
  
}
