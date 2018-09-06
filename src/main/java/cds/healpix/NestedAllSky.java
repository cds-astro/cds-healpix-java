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

import static cds.healpix.HealpixNestedBMOC.buildValue;
import static java.lang.Math.PI;

final class NestedAllSky implements HealpixNestedFixedRadiusConeComputer {

  private final double rRad;
  private final HealpixNestedBMOC moc;
  
  public NestedAllSky(final double radiusRad, final int depth) {
    assert radiusRad >= PI;
    this.rRad = radiusRad;
    final long[] mocVals = new long[12];
    for (int i = 0; i < 12; i++) {
      mocVals[i] = buildValue(0, i, true, depth);
    }
    this.moc = HealpixNestedBMOC.createUnsafe(depth, mocVals);
  }
  
  @Override
  public double getRadius() {
    return rRad;
  }

  @Override
  public HealpixNestedBMOC overlappingCells(double coneCenterLonRad, double coneCenterLatRad) {
    return this.moc;
  }

  @Override
  public HealpixNestedBMOC overlappingCenters(double coneCenterLonRad, double coneCenterLatRad) {
    return this.moc;
  }

  @Override
  public HealpixNestedFixedRadiusConeComputer newComputer() {
    return this;
  }

}
