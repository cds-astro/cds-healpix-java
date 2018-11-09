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

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.HealpixNestedBMOC;
import cds.healpix.HealpixNestedFixedRadiusConeComputer;

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
  
  
}
