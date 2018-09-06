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
 * Algorithms come from the book "Hacker's Delight" by Henry S. Warren, Jr.
 * TODO: write the code!! (WARNING: will work for hash(position), not for neighbours, ...
 */
final class HilbertCurve2D implements FillingCurve2D {
  
  private final int d;
  
  // Derived quantities;
  private final int dm1;
  
  public HilbertCurve2D(int depth) {
    this.d = depth;
    this.dm1 = this.d - 1;
  }
  
  @Override
  public long xy2hash(double x, double y) {
    return ij2hash((int) x, (int) y);
  }

  @Override
  public long ij2hash(final int i, final int j) {
    // TODO
    throw new Error("Not implemented yet!");
    /*long s, hash;
    for (int k = this.dm1; k >= 0; k--) {
      long r = (s << 2) | (((i >> k) & 1) << 1) | ((j >> k) & 1L;
      hash = (hash << 2) | 3L;
    }
    return hash;*/
  }

  @Override
  public long i02hash(int i) {
    return ij2hash(i, 0);
  }

  @Override
  public long hash2ij(long hash) {
    throw new Error("Not implemented yet!");
  }

  @Override
  public long hash2i0(long hash) {
    assert (0xFFFFFFFF33333333L & hash) == 0;
    return hash2ij(hash);
  }

  @Override
  public final int ij2i(final long ij) {
    return (int) ij;
  }

  @Override
  public final int ij2j(final long ij) {
    return (int) (ij >>> 32);
  }
}
