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

import cds.healpix.Healpix;
import cds.healpix.HealpixNestedFixedRadiusConeComputer;



public class NestedConeComputerApproxTest {

	
  @Test
  public void tesCone1() {
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
    // System.out.println("BMOC size: " + bmoc.size());
    
  }
	
}
