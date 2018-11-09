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

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.HealpixNested;
import cds.healpix.HealpixNestedBMOC;
import cds.healpix.NestedSmallCell;
import cds.healpix.RegularConeOrdinalHashComputer;

import static cds.healpix.Healpix.getBestStartingDepth;


public final class NestedSmallCellTest {

  
  public void coneTest(double coneCenterLonDeg, double coneCenterLatDeg, double radiusDeg, int depth,
      final long[] expectedRes) {
    /*System.out.println("Cone test. Center "
      + "lon: "  + coneCenterLonDeg + "deg, "
      + "lat: "  + coneCenterLatDeg + "deg, "
      + "radius: " + radiusDeg + " deg, "
      + "depth: " + depth
      );*/
    final double coneCenterLonRad = Math.toRadians(coneCenterLonDeg);
    final double coneCenterLatRad = Math.toRadians(coneCenterLatDeg);
    final double radiusRad = Math.toRadians(radiusDeg);
    final int startingDepth = getBestStartingDepth(radiusRad);
    final NestedSmallCell nlc = new NestedSmallCell(startingDepth, depth, radiusRad, RegularConeOrdinalHashComputer.UI);
    final HealpixNestedBMOC moc = nlc.overlappingCells(coneCenterLonRad, coneCenterLatRad);
    int i = 0;
    if (expectedRes == null) {
        System.out.println("MOC. size: " + moc.size() + "; deepSize: " + moc.computeDeepSize());
    } else {
      for (final HealpixNestedBMOC.CurrentValueAccessor cell : moc) {
        // System.out.println(cell);
        assertEquals(expectedRes[i++], cell.getHash());
      }
      assertEquals(expectedRes.length, i);
    }
  }
  
  
  // @Test
  public void latOfTest() {
    final int depth = 14;
    final HealpixNested hn = Healpix.getNested(depth);
    final double radiusRad = Math.toRadians(1 / 3600.0);
    final int startingDepth = getBestStartingDepth(radiusRad);
    final NestedSmallCell nsc = new NestedSmallCell(startingDepth, depth, radiusRad, RegularConeOrdinalHashComputer.UI);
    final int nRings = Healpix.nIsolatitudeRings(depth);
    assertEquals(0.5 * Math.PI, nsc.latOf(1, -1), 1e-16);
    assertEquals(-0.5 * Math.PI, nsc.latOf(1, nRings), 1e-16);
    /*System.out.println("iRing,latDeg");
    for (int i = 0; i < nRings; i++) {
      System.out.println(i + "," + Math.toDegrees(nsc.latOf(i)));
    }*/
  }

  
  @Test
  public void coneTestNorthPolarCap3() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(055.75724, +89.20594, 37.04')
//    coneTest(055.75724, +89.20594, 37.04 / 60d, 6, new long[]{4093, 4094, 4095, 8190, 8191});

//    coneTest(085.33262, -02.11126, 19.7 / 60d, 6, new long[]{22854, 22855});
    
//    coneTest(085.33262, -02.11126, 19.7 / 60d, 16, null);
  
    // draw phot(06:07:29.46,+41:36:38.7,2.018�)
    //coneTest(Math.toDegrees(1.6034818509470208), Math.toDegrees(0.7262444212606944), Math.toDegrees(0.03521504914227477), 13, null);
    // Cone. depth: 13; lonRad: 1.6034818509470208; latRad: 0.7262444212606944; rRad: 0.03521504914227477
    // Cone FX computed in 8.138862 ms
    // Regarder 13/111977291
    
    // Cone. depth: 11; lonRad: 3.5857899627457517; latRad: 1.5337734449159108; rRad: 0.10269986958660837
    // coneTest(Math.toDegrees(3.5857899627457517), Math.toDegrees(1.5337734449159108), Math.toDegrees(0.10269986958660837), 11, null);
    
    // draw phot(11:39:06.56,+89:55:53.8,16.21')
    // Cone. depth: 16; lonRad: 3.050439641567099; latRad: 1.5696028930332806; rRad: 0.004714625700186163
    coneTest(Math.toDegrees(3.050439641567099), Math.toDegrees(1.5696028930332806), Math.toDegrees(0.004714625700186163), 11, null);
    
    //Cone. depth: 9; lonRad: 1.9604995973603125; latRad: -1.4755999386429406; rRad: 0.43185443737294793
    coneTest(Math.toDegrees(1.9604995973603125), Math.toDegrees(-1.4755999386429406), Math.toDegrees(0.43185443737294793), 9, null);

    
    // draw circle(06:07:29.46,+41:36:38.7,2.018�)
    // draw phot(05:51:59.64,+84:06:52.9,1.317�)
  }
  
}
