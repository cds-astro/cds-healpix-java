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

package cds.healpix.common.math;

import org.junit.Test;
import static org.junit.Assert.*;

import cds.healpix.common.math.HackersDelight;

public final class HasckerDelightTestPerf {

  
  static void halfOfTest(int nElems) {
    final double[] data = new double[nElems];
    for (int i = 0; i < nElems; i++) {
      data[i] = java.lang.Math.random() * Math.PI;
    }
    long l1, l2;
    double s1, s2;
    
    l1 = System.nanoTime();
    s1 = halfOfBits(data);
    l2 = System.nanoTime();
    System.out.println("HalfByBits : " + s1 + " in " + (l2 - l1) / 1000000.0 + " ms");
    
    l1 = System.nanoTime();
    s2 = halfOfClassic(data);
    l2 = System.nanoTime();
    System.out.println("HalfClassic: " + s1 + " in " + (l2 - l1) / 1000000.0 + " ms");
    
    assertEquals(s1, s2);
    
    l1 = System.nanoTime();
    s1 = halfOfBits(data);
    l2 = System.nanoTime();
    System.out.println("HalfByBits : " + s1 + " in " + (l2 - l1) / 1000000.0 + " ms");
    
    l1 = System.nanoTime();
    s2 = halfOfClassic(data);
    l2 = System.nanoTime();
    System.out.println("HalfClassic: " + s1 + " in " + (l2 - l1) / 1000000.0 + " ms");
    
    assertEquals(s1, s2);
  }
  
  private static double halfOfClassic(final double[] values) {
    double s = 0;
    for (final double d : values) {
      s += 0.5 * d;
    }
    return s;
  }
  
  private static double halfOfBits(final double[] values) {
    double s = 0;
    for (final double d : values) {
      s += HackersDelight.halfOf(d);
    }
    return s;
  }
  
  @Test
  public void testPerf(final String[] args) {
    // Same exec time: smart compilo ;)
    System.out.println("Check that performances are the same (smart compilo ;) ):");
    halfOfTest(60000000);
  }
  
}
