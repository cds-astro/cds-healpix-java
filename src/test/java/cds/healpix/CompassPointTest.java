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

import static org.junit.Assert.*;

import org.junit.Test;

import cds.healpix.CompassPoint.MainWind;

public final class CompassPointTest {
  
  @Test
  public void mainWindsIndexTest() {
    // ^ SouthToWest axis
    // | _W(0,2) NW(1,2) _N(2,2)
    // | SW(0,1) _C(1,1) NE(2,1)
    // | _S(0,0) SE(1,0) _E(2,0)
    // -----------------------> SouthToEast axis
    
    assertEquals(MainWind.S , MainWind.getFromCoo(0, 0));
    assertEquals(MainWind.SE, MainWind.getFromCoo(1, 0));
    assertEquals(MainWind.E , MainWind.getFromCoo(2, 0));
    
    assertEquals(MainWind.SW, MainWind.getFromCoo(0, 1));
    assertEquals(MainWind.C , MainWind.getFromCoo(1, 1));
    assertEquals(MainWind.NE, MainWind.getFromCoo(2, 1));
    
    assertEquals(MainWind.W , MainWind.getFromCoo(0, 2));
    assertEquals(MainWind.NW, MainWind.getFromCoo(1, 2));
    assertEquals(MainWind.N , MainWind.getFromCoo(2, 2));
  }
  
}
