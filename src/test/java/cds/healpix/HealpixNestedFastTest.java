// Copyright 2017-2018 - Universite de Strasbourg/CNRS
// The CDS HEALPix library is developped by the Centre de Donnees
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

import static org.junit.Assert.*;

import org.junit.Test;

public final class HealpixNestedFastTest {

  
  @Test
  public void centerTest(){
    final int depth = 8;
    final long nCells = Healpix.nHash(depth);
    final HealpixNestedFast hnf = Healpix.getNestedFast(depth);
    final VerticesAndPathComputer vpc = hnf;
    
    VerticesAndPathComputer vpc2 = Healpix.getNested(depth);
    
    final HashComputer hc = hnf;
    final double[] lonlat = new double[2];
    final double[] lonlat2 = new double[2];
    long h2;
    for (long h = 0; h < nCells; h++) {
      // System.out.println("---- h = " + h + " ----");
      vpc.center(h, lonlat);
      vpc2.center(h, lonlat2);
      // System.out.println("lon: " + lonlat2[0] + "; lat: " + lonlat2[1]);
      h2 = hc.hash(lonlat[0], lonlat[1]);
      if (h != h2) {
        System.out.println("h2");
      }
      if (Math.abs(lonlat[0] - lonlat2[0]) > 1e-15) {
        System.out.println("lon is !=");
      }
      if (Math.abs(lonlat[1] - lonlat2[1]) > 1e-15) {
        System.out.println("lat is !=");
      }
      assertEquals(h, h2);
    }
  }
  
  @Test
  public void testCornerCase(){
    final int depth = 11;
    final double lon = 0.5 * Math.PI;
    final double lat = -1.1617;
    
    final HealpixNestedFast hnf = Healpix.getNestedFast(depth);
    System.out.println("ra dec:" + Math.toDegrees(lon) + " " + Math.toDegrees(lat));
    System.out.println(hnf.hash(lon, lat));
    
    final HealpixNested hn = Healpix.getNested(depth);
    System.out.println(hn.hash(lon, lat));
    
  }

}
