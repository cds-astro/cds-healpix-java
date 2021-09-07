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

  @Test
  public void tesMBT() {
    int order = 4;
    int hash = 87;
    CompassPoint.Cardinal cardinal = CompassPoint.Cardinal.N;
    VerticesAndPathComputer vc = Healpix.getNestedFast(order);
    double[] lonlat = vc.vertex(hash, cardinal);
    assert true;
  }

  @Test
  public void tesMBT2() {
    int order = 1;
    int hash = 3;
    CompassPoint.Cardinal cardinal = CompassPoint.Cardinal.N;
    VerticesAndPathComputer vc = Healpix.getNestedFast(order);
    double[] lonlat = vc.vertex(hash, cardinal);
    assert Math.toDegrees(lonlat[0]) == 45.0;
    assert Math.toDegrees(lonlat[1]) == 90.0;
  }

  @Test
  public void testAllVertices() {
    // Test that all asserts are ok and values are in expected bounds
    for(int o = 1; o < 10; o++) {
      VerticesAndPathComputer vc = Healpix.getNestedFast(o);
      for(long h = 0; h < 12 * (1L << (o << 1)) ; h++) {
        double[] lonlatN = vc.vertex(h,  CompassPoint.Cardinal.N);
        checklonlat(o, h, 'N', lonlatN);
        double[] lonlatS = vc.vertex(h,  CompassPoint.Cardinal.S);
        checklonlat(o, h, 'S', lonlatS);
        double[] lonlatE = vc.vertex(h,  CompassPoint.Cardinal.E);
        checklonlat(o, h, 'E', lonlatE);
        double[] lonlatW = vc.vertex(h,  CompassPoint.Cardinal.W);
        checklonlat(o, h, 'W', lonlatW);
      }
    }
  }

  void checklonlat(int o, long h, char c, double[] lonlat) {
    assert lonlat[0] >= 0 && lonlat[0] <= 2.0 * Math.PI
        : "o: " + o + "; h: " + h + "; dir: " + c + "; lon:" + lonlat[0];
    assert lonlat[1] >= -0.5 * Math.PI && lonlat[1] <= 0.5 * Math.PI
        : "o: " + o + "; h: " + h + "; dir: " + c + "; lat:" + lonlat[1];
  }

  @Test
  public void crossCheckVerticesAndVerticesFast() {
    for(int o = 1; o < 9; o++) {
      VerticesAndPathComputer vcr = Healpix.getNested(o);
      VerticesAndPathComputer vcf = Healpix.getNestedFast(o);
      for(long h = 0; h < 12 * (1L << (o << 1)) ; h++) {
        double[] lonlatN_regu = vcr.vertex(h,  CompassPoint.Cardinal.N);
        double[] lonlatN_fast = vcf.vertex(h,  CompassPoint.Cardinal.N);
        checkEq(o, h, 'N', lonlatN_regu, lonlatN_fast);
        double[] lonlatS_regu = vcr.vertex(h,  CompassPoint.Cardinal.S);
        double[] lonlatS_fast = vcf.vertex(h,  CompassPoint.Cardinal.S);
        checkEq(o, h, 'S', lonlatS_regu, lonlatS_fast);
        double[] lonlatE_regu = vcr.vertex(h,  CompassPoint.Cardinal.E);
        double[] lonlatE_fast = vcf.vertex(h,  CompassPoint.Cardinal.E);
        checkEq(o, h, 'E', lonlatE_regu, lonlatE_fast);
        double[] lonlatW_regu = vcr.vertex(h,  CompassPoint.Cardinal.W);
        double[] lonlatW_fast = vcf.vertex(h,  CompassPoint.Cardinal.W);
        checkEq(o, h, 'W', lonlatW_regu, lonlatW_fast);
      }
    }
  }

  void checkEq(int o, long h, char c, double[] lonlat1, double[] lonlat2) {
    assert Math.abs(lonlat1[0] - lonlat2[0]) < 1e-13
        || Math.abs(lonlat1[0] - lonlat2[0]) - 2 * Math.PI < 1e-13 // 2*PI possible because of numerical precision
      : "o: " + o + "; h: " + h + "; dir: " + c + "; lon: " + lonlat1[0] + " != " + lonlat2[0];
    assert Math.abs(lonlat1[1] - lonlat2[1]) < 1e-13
      : "o: " + o + "; h: " + h + "; dir: " + c + "; lat: " + lonlat1[1] + " != " + lonlat2[1];
  }
  
}
