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
import static org.junit.Assert.assertTrue;

import java.util.logging.Level;

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.HealpixNestedBMOC;
import cds.healpix.HealpixNestedFixedRadiusConeComputer;
import cds.healpix.NestedSmallCellApproxedMethod.GrowableLongArray;

public class NestedSmallCellApproxedMethodTest {


  public void coneTest(double coneCenterLonRad, double coneCenterLatRad, double radiusRad, int depth,
      final long[] expectedRes) {
    // final int startingDepth = Constants.getBestStartingDepthForConeSearch(radiusRad);
    // final NestedLargeCellApproxedMethod nlc = new NestedLargeCellApproxedMethod(startingDepth, depth, radiusRad);
    // final NestedSmallCellApproxedMethod nlc = new NestedSmallCellApproxedMethod(startingDepth, depth, radiusRad);
    final HealpixNestedFixedRadiusConeComputer nlc = Healpix.getNested(depth).newConeComputerApprox(radiusRad);

    System.out.println(nlc.getClass().getName());

    final HealpixNestedBMOC moc = nlc.overlappingCells(coneCenterLonRad, coneCenterLatRad);
    System.out.println("Moc size: " + moc.size());
    System.out.println("Moc deep size: " + moc.computeDeepSize());
    int i = 0;
    for (final HealpixNestedBMOC.CurrentValueAccessor cell : moc) {
      // System.out.println(cell);
      // assertEquals(expectedRes[i++], cell.getHash());
    }
    // assertEquals(expectedRes.length, i);
  }

  @Test
  public void testCone() {
    System.out.println("lon: " + Math.toRadians(085.51340) + "; lat: " + Math.toRadians(-01.89284));
    //085.51340 -01.89284

    //coneTest(0.0, 0.0, 2.0, 3, null);
    // Cone. depth: 6; lonRad: 1.4887397926478119; latRad: -0.03928177631283567; rRad: 0.01994913445007241
    // coneTest(1.4887397926478119, -0.03928177631283567, 0.01994913445007241, 6, null);
    System.out.println("---------");
    coneTest(1.4923946835544672, -0.03320059031799468, 0.016729491874516444, 7, null);
    // Cone. depth: 7; lonRad: 1.4923946835544672; latRad: -0.03320059031799468; rRad: 0.016729491874516444

    // Cone. depth: 3; lonRad: 1.4985159126238619; latRad: 1.4771883195119886; rRad: 0.642057147736403
    coneTest(1.4985159126238619, 1.4771883195119886, 0.642057147736403, 3, null);
  }



  @Test
  public void tesConeSeeGithubIssue15() {
    // See https://github.com/cds-astro/cds-healpix-java/issues/15
    int depth = 14;
    double radius = 0.001;
    double alpha = 0.002;
    double delta = -1.3;
    HealpixNestedFixedRadiusConeComputer coner =
        Healpix.getNested(depth).newConeComputerApprox(radius);
    final HealpixNestedBMOC bmoc = coner.overlappingCells(alpha, delta);
    // Test was failling (out of bound exception) due to a too small initial 
    // MOC size guess.
    // Not assert, ok if does not fail any more.    
  }

  // Code from Mark Taylor, see 
  // https://github.com/cds-astro/cds-healpix-java/issues/15
  private boolean subTest(int capacity, int n) {
    NestedSmallCellApproxedMethod.GrowableLongArray list = new GrowableLongArray(capacity);
    for (int i = 0; i < n; i++) {
      list.add(i);
    }
    long[] array = list.getArray();
    int count = list.getCursor();

    boolean ok = true;
    ok = ok && (count == n);
    for (int i = 0; i < n; i++) {
      ok = ok && ( array[i] == i );
    }
    /*if (!ok) {
          throw new RuntimeException("Failure");
      }*/
    return ok;
  }

  // Code from Mark Taylor, see 
  // https://github.com/cds-astro/cds-healpix-java/issues/15
  @Test
  public void growableLongArrayTest() {
    // Remove spurious message during these specific tests
    NestedSmallCellApproxedMethod.GrowableLongArray.LOGGER.setLevel(Level.SEVERE);
    for (int ic = 10; ic < 100; ic++) {
      for (int in = 0; in < 100; in++) {
        assertTrue(this.subTest(ic, in));
      }
    }
    NestedSmallCellApproxedMethod.GrowableLongArray.LOGGER.setLevel(Level.WARNING);
  }


}
