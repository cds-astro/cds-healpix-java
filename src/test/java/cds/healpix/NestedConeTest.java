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

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.HealpixNestedBMOC;
import cds.healpix.HealpixNestedFixedRadiusConeComputer;

public class NestedConeTest {

  public void coneTest(double coneCenterLonDeg, double coneCenterLatDeg, double radiusDeg, int depth,
      final long[] expectedRes) {
    System.out.println("Cone test. Center "
        + "lon: "  + coneCenterLonDeg + "deg, "
        + "lat: "  + coneCenterLatDeg + "deg, "
        + "radius: " + radiusDeg + " deg, "
        + "depth: " + depth
        );
    final double coneCenterLonRad = Math.toRadians(coneCenterLonDeg);
    final double coneCenterLatRad = Math.toRadians(coneCenterLatDeg);
    final double radiusRad = Math.toRadians(radiusDeg);
    final HealpixNestedFixedRadiusConeComputer cc = Healpix.getNested(depth).newConeComputer(radiusRad);
    final HealpixNestedBMOC moc = cc.overlappingCells(coneCenterLonRad, coneCenterLatRad);
    int i = 0;
    if (expectedRes != null) {
      for (final HealpixNestedBMOC.CurrentValueAccessor cell : moc) {
        // System.out.println(cell);
        assertEquals(expectedRes[i++], cell.getHash());
      }
      assertEquals(expectedRes.length, i);
    } else {
      for (final HealpixNestedBMOC.CurrentValueAccessor cell : moc) {
        System.out.println(cell);
      }
    }
  }
  
  public void coneCentersTest(double coneCenterLonDeg, double coneCenterLatDeg, double radiusDeg, int depth,
      final long[] expectedRes) {
    System.out.println("Cone test. Center "
        + "lon: "  + coneCenterLonDeg + "deg, "
        + "lat: "  + coneCenterLatDeg + "deg, "
        + "radius: " + radiusDeg + " deg, "
        + "depth: " + depth
        );
    final double coneCenterLonRad = Math.toRadians(coneCenterLonDeg);
    final double coneCenterLatRad = Math.toRadians(coneCenterLatDeg);
    final double radiusRad = Math.toRadians(radiusDeg);
    final HealpixNestedFixedRadiusConeComputer cc = Healpix.getNested(depth).newConeComputer(radiusRad);
    final HealpixNestedBMOC moc = cc.overlappingCenters(coneCenterLonRad, coneCenterLatRad);
    int i = 0;
    if (expectedRes != null) {
      for (final HealpixNestedBMOC.CurrentValueAccessor cell : moc) {
        // System.out.println(cell);
        assertEquals(expectedRes[i++], cell.getHash());
      }
      assertEquals(expectedRes.length, i);
    } else {
      for (final HealpixNestedBMOC.CurrentValueAccessor cell : moc) {
        // System.out.println(cell);
        // if (cell.getHash() == 373286564231L) System.out.println("COUCOU");
      }
    }
  }

  double lonFromSexa(final String lon) {
    final String[] elems = lon.trim().split(" +");
    final double hour = Integer.parseInt(elems[0]);
    final double minh = Integer.parseInt(elems[1]) / 60.0;
    final double sech = Double.parseDouble(elems[2]) / 3600.0;
    return 15 * (hour + minh + sech);
  }

  double latFromSexa(final String lat) {
    final String[] elems = lat.trim().split(" +");
    final double deg = Integer.parseInt(elems[0]);
    final double min = Integer.parseInt(elems[1]) / 60.0;
    final double sec = Double.parseDouble(elems[2]) / 3600.0;
    return deg > 0 ? deg + min + sec : deg - min - sec; 
  }

  @Test
  public void tesConeNP() {
    // draw circle(61.0749875, +89.90546338888889, 5.2403223034772175 deg)
    double lon = lonFromSexa("04 04 17.99700");
    double lat = latFromSexa("+89 54 19.6682");
    double rad = 5.2403223034772175;
    int order = 3;
    coneTest(lon, lat, rad, order, new long[]{63, 127, 191, 255});
  }
  
  @Test
  public void tesConeNP2() {
    // draw circle(37.954541666666664, +89.2641111111111, 0.8410004718589024 deg)
    double lon = lonFromSexa("02 31 49.09000");
    double lat = latFromSexa("+89 15 50.8000");
    double rad = 0.8410004718589024;
    int order = 6;
    coneTest(lon, lat, rad, order, null); //new long[]{4092, 4093, 4094, 4095, 8190, 8191, 12287, 16381, 16383});
  }
  
  @Test
  public void tesCone3() {
    // draw circle(37.954541666666664, +89.2641111111111, 0.8410004718589024 deg)
    double lon = lonFromSexa("06 30 12.40");
    double lat = latFromSexa("+04 50 26.7");
    double rad = 30.5 / 3600d;
    int order = 18;
    coneCentersTest(lon, lat, rad, order, null); //new long[]{4092, 4093, 4094, 4095, 8190, 8191, 12287, 16381, 16383});
    // 18/373286564231 dedans?
  }

  public static void main(String[] args) {
    new NestedConeTest().tesCone3();
  }
  
}
