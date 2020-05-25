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

import cds.healpix.HashComputer;
import cds.healpix.Healpix;
import cds.healpix.HealpixNested;

public class HealpixNestedHashComputerTest {

  @Test
  public void test() {
    HealpixNested hn = Healpix.getNested(0);
    HashComputer hc = hn.newHashComputer();
    
    long h1 = hc.hash(Math.toRadians(85.30250995302792), Math.toRadians(-2.2546799530640813));
    long h2 = hc.hash(1.4888095533470505, -0.03935148329243798);
    long h3 = hc.hash(1.4888095533470505,  0.602705664443962);
    long h4 = hc.hash(1.4888095533470505, -0.681408631028838);
    long h5 = hc.hash(1.4888095888031632, -0.03935151872124723);
    long h6 = hc.hash(0.6641163101020234, 1.559806382359946);
    long h7 = hc.hash(3.8961232510963066, 1.5707963267948966);
    long h8 = hc.hash(0.7545305975065133, 1.5707963267948966);
    long h9 = hc.hash(5.955948491951194, 1.5528603751397578);
    long h10 =  hc.hash(6.841839252888341, 1.5707963267948966);
    long h11 = hc.hash(1.4985159126238619, 1.4771883195119886);

    assertEquals(5L, h1);
    assertEquals(5L, h2);
    assertEquals(5L, h3);
    assertEquals(8L, h4);
    assertEquals(5L, h5);
    assertEquals(0L, h6);
    assertEquals(2L, h7);
    assertEquals(0L, h8);
    assertEquals(3L, h9);
    assertEquals(0L, h10);
    assertEquals(0L, h11);
  }
  
  
  
}
