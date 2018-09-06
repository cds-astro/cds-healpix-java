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
 * 
 * In this class, the projetion ranges are:
 * (lon, lat) ----> (x, y)
 * n in [0, 10]
 * lon in [-n * PI, n * PI]
 * lat in [-PI/2, PI/2]
 * x in [-4 * n, 4 * n]
 * y in [-2, 2]
 * 
 * @author F.-X. Pineau
 *
 */
class HealpixProjection implements Projection {

  private final HealpixProjector projector = new HealpixProjector();
  private final HealpixUnprojector unprojector = new HealpixUnprojector();
  
  @Override
  public double[] project(final double lon, final double lat) {
    final double[] xy = new double[2];
    project(lon, lat, xy);
    return xy;
  }

  @Override
  public void project(final double lon, final double lat, final double[] resultXY) {
    projector.project(lon, lat, resultXY);
  }
  
  @Override
  public double[] unproject(final double x, final double y) {
    final double[] lonlat = new double[2];
    unproject(x, y, lonlat);
    return lonlat;
  }

  @Override
  public void unproject(final double x, final double y, final double[] resultLonLat) {
    unprojector.unproject(x, y, resultLonLat);
  }
}
