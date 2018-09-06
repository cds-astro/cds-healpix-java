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

import cds.healpix.FlatHashIterator;
import cds.healpix.HealpixNestedBMOC;
import cds.healpix.HealpixNestedBMOC.CurrentValueAccessor;

public class HealpixBMOCTest {

  @Test
  public void createPackingTest1() {
    final long[] rawMoc = new long[8];
    rawMoc[0] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0001" + "10" + "10", 2), true, 3);
    rawMoc[1] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0001" + "10" + "11", 2), false, 3);
    
    rawMoc[2] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "00", 2), true, 3);
    rawMoc[3] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "01", 2), true, 3);
    rawMoc[4] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "10", 2), true, 3);
    rawMoc[5] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "11", 2), true, 3);
    
    rawMoc[6] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0101" + "10" + "00", 2), true, 3);
    rawMoc[7] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0111" + "10" + "10", 2), false, 3);
    
    HealpixNestedBMOC moc = HealpixNestedBMOC.createPacking(3, rawMoc, 8);

    int i = 0;
    assertEquals(5, moc.size());
    for (final CurrentValueAccessor cva : moc) {
      // System.out.println("depth: " + cva.getDepth() + "; val: " + String.format("%8s", Long.toBinaryString(cva.getHash())) + "; full: " + cva.isFull());
      switch(i) {
      case 0: assertEquals(rawMoc[0], cva.getRawValue()); break;
      case 1: assertEquals(rawMoc[1], cva.getRawValue()); break;
      case 2: assertEquals(HealpixNestedBMOC.buildValue(1, Long.parseLong("0011" + "10", 2), true, 3), cva.getRawValue()); break;
      case 3: assertEquals(rawMoc[6], cva.getRawValue()); break;
      case 4: assertEquals(rawMoc[7], cva.getRawValue()); break;
      default: assertEquals(false, true);
      }
      i++;
    }
    assertEquals(32, moc.computeDeepSize());
    
    FlatHashIterator fit = moc.flatHashIterator();
    i = 0;
    while (fit.hasNext()) {
       fit.next();
       i++;
    }
    assertEquals(32, i);
  }
  
  @Test
  public void createPackingTest2() {
    final long[] rawMoc = new long[20];
    rawMoc[0] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0001" + "10" + "10", 2), true, 3);
    rawMoc[1] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0001" + "10" + "11", 2), false, 3);
    
    rawMoc[2] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "00" + "00", 2), true, 3);
    rawMoc[3] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "00" + "01", 2), true, 3);
    rawMoc[4] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "00" + "10", 2), true, 3);
    rawMoc[5] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "00" + "11", 2), true, 3);
    
    rawMoc[6] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "01" + "00", 2), true, 3);
    rawMoc[7] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "01" + "01", 2), true, 3);
    rawMoc[8] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "01" + "10", 2), true, 3);
    rawMoc[9] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "01" + "11", 2), true, 3);
    
    rawMoc[10] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "00", 2), true, 3);
    rawMoc[11] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "01", 2), true, 3);
    rawMoc[12] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "10", 2), true, 3);
    rawMoc[13] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "10" + "11", 2), true, 3);
    
    rawMoc[14] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "11" + "00", 2), true, 3);
    rawMoc[15] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "11" + "01", 2), true, 3);
    rawMoc[16] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "11" + "10", 2), true, 3);
    rawMoc[17] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0011" + "11" + "11", 2), true, 3);
    
    rawMoc[18] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0101" + "10" + "00", 2), true, 3);
    rawMoc[19] = HealpixNestedBMOC.buildValue(2, Long.parseLong("0111" + "10" + "10", 2), false, 3);
    
    HealpixNestedBMOC moc = HealpixNestedBMOC.createPacking(3, rawMoc);

    int i = 0;
    assertEquals(5, moc.size());
    for (final CurrentValueAccessor cva : moc) {
      switch(i) {
      case 0: assertEquals(rawMoc[0], cva.getRawValue()); break;
      case 1: assertEquals(rawMoc[1], cva.getRawValue()); break;
      case 2: assertEquals(HealpixNestedBMOC.buildValue(0, Long.parseLong("0011", 2), true, 3), cva.getRawValue()); break;
      case 3: assertEquals(rawMoc[18], cva.getRawValue()); break;
      case 4: assertEquals(rawMoc[19], cva.getRawValue()); break;
      default: assertEquals(false, true);
      }
      i++;
    }
  }
  
}
