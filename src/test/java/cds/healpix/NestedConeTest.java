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
    System.out.println(cc.getClass().getName());
    
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
        System.out.println(cell);
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
    coneTest(lon, lat, rad, order, new long[]{4092, 4093, 4094, 4095, 8190, 8191, 12287, 16381, 16383});
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

  @Test
  public void testCone4() {
    // draw circle(0.0, +0.0, 37.0 deg)
    double lon = 0;
    double lat = 0;
    double rad = 37;
    int order = 3;
   
    coneTest(lon, lat, rad, order, null);
    
    //final HealpixNested hn = Healpix.getNested(order);
    //final HealpixNestedFixedRadiusConeComputer cp = hn.newConeComputer( Math.toRadians(radius) );
    //final HealpixNestedBMOC bmoc = cp.overlappingCells(Math.toRadians(ra), Math.toRadians(dec));
    // long [] out = toFlatArrayOfHash(bmoc);
   
    // System.out.print("overlappingCells ra="+ra+" dec="+dec+" radius="+radius+" order="+order+"\n=> ");
    // for( long a : out ) System.out.print(" "+a);
    // System.out.println();

  }
  
  
  @Test
  public void testCone5() {
    double lon = 012.46214169;
    double lat = -72.49285924;
    double rad = 2.6 * 1.42;
    int order = 4;
    // draw circle(12.46214169,+72.49285924,3.6919999999999997)
    // draw moc 4/2058 2059 2080 2081 2082 2083
    coneTest(lon, -lat, rad, order, new long[]{233, 235, 236, 237, 238, 239, 248, 250});
    // draw circle(12.46214169,-72.49285924,3.6919999999999997)
    // draw moc 4/2058 2059 2080 2081 2082 2083
    coneTest(lon, lat, rad, order, new long[]{2058, 2059, 2080, 2081, 2082, 2083, 2088, 2089});
   
  }

  @Test
  public void testCone6() {
    // double lon = lonFromSexa("00 55 58.86840");
    // double lat = latFromSexa("-73 22 13.6829");
    double lon = lonFromSexa("00 53 35.72");
    double lat = latFromSexa("-72 57 44.8");
    double rad = 4.6;
    int order = 4;
    coneTest(lon, lat, rad, order, new long[]{2058, 2059, 2080, 2081, 2082, 2083, 2088, 2089, 2832, 2833, 2836});
    // draw circle( 00:55:58.86840, -73:22:13.6829, 4.6)
    // draw moc 4/2058 2059 2080 2081 2082 2083 2088 2089 2832 2833 2836*/
   
  }

  
  @Test
  public void testCone7() {
    double lon = 13.158329;
    double lat = -72.80028;
    double rad = 5.64323;
    int order = 3;
    
    // draw circle(13.158329,+72.80028,5.64323)
    // draw circle(103.158329,+72.80028,5.64323)
    coneTest(lon, -lat, rad, order, new long[]{57, 58, 59, 60, 62, 245, 247});
    
    // draw circle(13.158329,-72.80028,5.64323)
    // draw moc 3/ 514 515 520 521 522
    coneTest(lon, lat, rad, order, new long[]{514, 515, 520, 521, 522, 708, 709});

  }
  
  @Test
  public void testCone8() {
    double lon = 97.91750000;
    double lat = 5.76952778;
    double rad = 4.45 / 60;
    int order = 9;
   
    // draw circle(97.91750000,+5.76952778, 0.07416666666666666666)
    coneTest(lon, lat, rad, order, new long[]{1424261, 1424262, 1424263, 1424269, 1424274, 1424280});
    
    /*final HealpixNested hn = Healpix.getNested(order);
    final HealpixNestedFixedRadiusConeComputer cp = hn.newConeComputer( Math.toRadians(radius) );
    final HealpixNestedBMOC bmoc = cp.overlappingCells(Math.toRadians(ra), Math.toRadians(dec));
    long [] out = toFlatArrayOfHash(bmoc);
   
    System.out.print("overlappingCells checker:\ndraw circle("+ra+","+dec+","+radius+")\ndraw moc "+order+"/");
    for( long a : out ) System.out.print(" "+a);
    System.out.println();*/
  }
  
  public static void main(String[] args) {
    // new NestedConeTest().tesCone3();
    new NestedConeTest().testCone8();
  }
  
}
