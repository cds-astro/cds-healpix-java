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
import cds.healpix.HealpixNested;
import cds.healpix.HealpixNestedBMOC.CurrentValueAccessor;
import cds.healpix.HealpixNestedPolygonComputer;
import cds.healpix.common.math.Math;

public final class HealpixNestedPolygonComputerTest {

  
  // draw polygon(10.0, 0.0, 0.0, -2.0, 2.0, 10.0, 10.0, 0.0, 2.5, 0.0, 4.0, 5.0, 10.0, 0.0, 3.0, 0.0 ,3.5, 1.0, 4.5, 0.0) 
  
  // ERREUR:
// draw polygon(03:15:19.13,-61:24:30.1, 00:20:42.90,-68:26:23.1, 23:03:34.78,-57:10:39.9, 23:03:34.78,-57:10:39.9, 08:49:10.85,-18:22:49.4, 21:09:31.59,-47:21:14.9, 21:06:03.04,-16:23:58.5, 20:31:13.83,-01:20:21.4, 21:06:19.28,+18:52:58.0, 21:06:19.28,+18:52:58.0, 21:19:38.62,+21:35:25.7, 21:25:33.81,+22:52:13.3, 22:17:12.93,+73:04:48.1, 22:17:12.93,+73:04:48.1, 22:09:01.04,+75:26:21.5, 13:58:14.81,+78:06:11.9, 04:16:45.73,+00:23:33.3, 13:21:07.94,+56:39:05.6, 16:13:10.95,+18:23:40.9, 10:18:25.58,+20:16:58.0, 13:04:32.91,-10:44:20.5, 12:06:34.45,+30:46:57.7, 08:59:25.66,+11:01:33.1, 11:34:10.09,-12:59:23.6, 13:34:25.91,-26:55:19.6, 15:38:34.12,+05:11:07.6, 16:20:54.75,-32:39:02.0, 15:33:20.32,-63:33:48.7, 08:13:39.81,-68:19:53.3, 09:42:20.97,-42:05:28.2, 04:54:53.82,-60:33:08.0)
  
  /*draw polygon(
      05:41:16.32,-02:05:06.6, 
      05:42:01.24,-02:19:37.8,
      05:41:02.15,-02:27:25.8,
      05:40:21.32,-02:17:24.6,
      05:40:31.41,-02:06:15.0)*/
  
 // draw polygon(05:41:16.32,-02:05:06.6, 05:42:01.24,-02:19:37.8, 05:41:02.15,-02:27:25.8, 05:40:21.32,-02:17:24.6, 05:40:31.41,-02:06:15.0)
  
 // draw polygon(12:09:52.27,+88:48:06.3, 08:07:01.96,+88:50:31.4, 06:46:07.27,+89:09:22.9, 05:20:08.95,+88:26:11.9, 03:47:41.42,+88:39:33.6, 02:03:24.29,+87:49:41.7, 01:24:27.20,+88:24:58.6, 23:21:44.90,+88:23:06.1, 20:28:09.52,+89:07:08.7, 21:35:39.19,+89:18:26.3, 17:19:14.70,+89:09:47.0, 15:17:45.99,+89:08:02.9, 14:17:15.54,+88:28:27.2, 13:53:17.54,+89:00:26.5, 12:53:00.72,+88:22:56.4)
  
 // draw polygon(14:22:24.12,+14:09:33.6, 13:20:19.08,+02:01:08.1, 11:49:48.33,-09:43:32.6, 10:47:20.93,-14:52:07.9, 08:38:06.57,-22:06:58.9, 06:47:36.01,-28:17:21.6, 04:46:52.23,-28:35:23.0, 02:29:06.86,-28:41:31.0, 02:29:06.86,-28:41:31.0, 02:21:35.03,-28:35:48.2, 00:47:33.31,-27:09:20.3, 23:07:46.45,-23:37:28.5, 21:05:45.41,-16:48:35.8, 20:25:57.20,-07:49:42.6, 20:20:40.15,+08:08:01.7, 20:58:04.57,+09:53:53.4, 22:07:00.01,+02:23:27.2, 00:50:43.94,-03:13:59.3, 03:07:48.93,-04:19:24.9, 04:53:03.08,-02:07:36.7, 04:53:03.08,-02:07:36.7, 04:59:46.12,-01:25:20.7, 06:41:11.36,+04:06:05.7, 06:41:11.36,+04:06:05.7, 08:43:58.98,+11:21:41.1, 08:43:58.98,+11:21:41.1, 10:34:29.93,+18:55:07.2, 10:34:29.93,+18:55:07.2, 10:43:50.34,+19:37:25.6, 13:00:45.39,+23:04:32.7, 13:00:45.39,+23:04:32.7, 14:22:00.57,+23:01:24.1, 14:41:59.11,+18:57:44.5)

 /* @Test
  public void polygonTest1() {
    final String[] verticesSexa = new String[]{
        "12:09:52.27,+88:48:06.3", 
        "08:07:01.96,+88:50:31.4",
        "06:46:07.27,+89:09:22.9",
        "05:20:08.95,+88:26:11.9",
        "03:47:41.42,+88:39:33.6",
        "02:03:24.29,+87:49:41.7",
        "01:24:27.20,+88:24:58.6",
        "23:21:44.90,+88:23:06.1",
        "20:28:09.52,+89:07:08.7",
        "21:35:39.19,+89:18:26.3",
        "17:19:14.70,+89:09:47.0",
        "15:17:45.99,+89:08:02.9",
        "14:17:15.54,+88:28:27.2",
        "13:53:17.54,+89:00:26.5",
        "12:53:00.72,+88:22:56.4"
    };
    HealpixNested hn = Healpix.getNested(13);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    final double[][] vertices = new double[verticesSexa.length][2];
    for (int i = 0; i < verticesSexa.length; i++) {
      final String[] lonlatSexa = verticesSexa[i].split(",");
      final String[] lonSexa = lonlatSexa[0].split(":");
      final String[] latSexa = lonlatSexa[1].split(":");
      final double lon = 15 * (Integer.parseInt(lonSexa[0]) + Integer.parseInt(lonSexa[1]) / 60d + Double.parseDouble(lonSexa[2]) / 3600d);
      double lat = Integer.parseInt(latSexa[0]);
      if (lat < 0) {
        lat -= Integer.parseInt(latSexa[1]) / 60d + Double.parseDouble(latSexa[2]) / 3600d;
      } else {
        lat += Integer.parseInt(latSexa[1]) / 60d + Double.parseDouble(latSexa[2]) / 3600d;;
      }
      vertices[i][0] = Math.toRadians(lon);
      vertices[i][1] = Math.toRadians(lat);
    }
    pc.overlappingCells(vertices);
    
  }
  
  @Test
  public void polygonTest2() {
    final String[] verticesDec = new String[]{
        "287.18290328,-02.26765169",
        "287.19157232,-02.27546342",
        "287.18274564,-02.27946381",
        "287.18252498,-02.27562093",
        "287.17477012,-02.27549490"
    };
    HealpixNested hn = Healpix.getNested(18);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    final double[][] vertices = new double[verticesDec.length][2];
    for (int i = 0; i < verticesDec.length; i++) {
      final String[] lonlatDec = verticesDec[i].split(",");
      vertices[i][0] = Math.toRadians(Double.parseDouble(lonlatDec[0]));
      vertices[i][1] = Math.toRadians(Double.parseDouble(lonlatDec[1]));
    }
    
    final HealpixNestedBMOC bmoc = pc.overlappingCenters(vertices);
    
    System.out.println("MOC VIEW");
    for (CurrentValueAccessor curr : bmoc) {
      System.out.println(curr);
    }
    
    System.out.println("FLAT VIEW");
    final FlatHashIterator it = bmoc.flatHashIterator();
    while (it.hasNext()) {
      System.out.println(it.next());
    }
  }*/
  
  
  @Test
  public void polyonTest3() {
    final int depth = 3;
    final String[] verticesDec = new String[]{
        "0.0, 0.0",
        "0.0, 0.5",
        "0.25, 0.25"
    };
    final long[] expected = new long[]{304, 305, 306, 307, 308, 310, 313, 316};

    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    final double[][] vertices = new double[verticesDec.length][2];
    for (int i = 0; i < verticesDec.length; i++) {
      final String[] lonlatDec = verticesDec[i].split(",");
      vertices[i][0] = Double.parseDouble(lonlatDec[0]);
      vertices[i][1] = Double.parseDouble(lonlatDec[1]);
    }
    
    final HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      assertEquals(expected[i], it.next());
    }
    
    /*System.out.println("MOC VIEW:");
    for (CurrentValueAccessor curr : bmoc) {
      System.out.println(curr);
    }
    System.out.println("FLAT VIEW");
    final FlatHashIterator it = bmoc.flatHashIterator();
    while (it.hasNext()) {
      System.out.println(it.next());
    }*/
  }
  
  
  @Test
  public void polyonTest4() {
    // Test special points
    // Aladin: draw polygon(65.11781779000003, 85.012424, 89.70533626000001, 87.06130188, 60.23667431000001, 85.609882)
    // - input data
    final int depth = 6;
    final double[][] vertices = deg2rad(new double[][]{
      {65.11781779000003, 85.012424}, 
      {89.70533626000001, 87.06130188}, 
      {60.23667431000001, 85.609882}
    });
    final long[] expected = new long[]{4062, 4063, 4084, 4085};
    // - generic algo
    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      assertEquals(expected[i], it.next());
    }
  }
  
  @Test
  public void polyonTest5() {
    // Aladin: draw polygon(359.70533626,+87.06130188, 330.23667431,+85.60988200, 335.11781779,+85.01242400)
    // Test special points
    // - input data
    final int depth = 6;
    final double[][] vertices = deg2rad(new double[][]{
      {359.70533626, 87.06130188}, 
      {330.23667431, 85.60988200},
      {335.11781779, 85.01242400}
    });
    final long[] expected = new long[]{16350, 16351, 16372, 16373};
    // - generic algo
    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      assertEquals(expected[i], it.next());
    }
  }
  
  @Test
  public void polyonTest6() {
    // Aladin:  draw polygon(224.86211710,+78.10924662, 176.91129363 +83.92878811, 135.81578643,+78.24840426, 200.73574863,+73.58038790)
    // Test special points
    // - input data
    final int depth = 3;
    final double[][] vertices = deg2rad(new double[][]{
      {224.86211710, 78.10924662},
      {176.91129363, 83.92878811},
      {135.81578643, 78.24840426},
      {200.73574863, 73.58038790}
    });
    final long[] expected = new long[]{119, 125, 127, 187, 188, 190, 191};
    // - generic algo
    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      assertEquals(expected[i], it.next());
    }
  }
  
  @Test
  public void polyonTest7() {
    // Aladin: draw polygon(359.70533626,-87.06130188, 330.23667431,-85.60988200, 335.11781779,-85.01242400)
    // Test special points
    // - input data
    final int depth = 6;
    final double[][] vertices = deg2rad(new double[][]{
      {359.70533626, -87.06130188}, 
      {330.23667431, -85.60988200}, 
      {335.11781779, -85.01242400}
    });
    final long[] expected = new long[]{45061, 45063, 45072, 45074};
    // - generic algo
    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      assertEquals(expected[i], it.next());
    }
  }
  
  @Test
  public void polyonTest8() {
    // In Aladin: draw polygon(180.08758393,-41.43289179, 191.00310758,-29.99207687, 181.59160475,-34.21976170)
    // Test special points
    // - input data
    final int depth = 3;
    final double[][] vertices = deg2rad(new double[][]{
      {180.08758393,-41.43289179},
      {191.00310758,-29.99207687},
      {181.59160475,-34.219761700}
    });
    final long[] expected = new long[]{384, 385, 682, 683};
    // - generic algo
    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      assertEquals(expected[i], it.next());
    }
  }
  
  @Test
  public void polyonTest9() {
    // Aladin:
    /// draw polygon(0.8537903597029015,-49.60560772150546,359.7264963624239,-49.60560826882517,359.73477207820594,-48.875097433186816,0.8455129889279733,-48.875096899797896)
    // final int depth = 5; // final long[] expected = new long[]{8840, 8842, 11588, 11589};
    final int depth = 10;
    /*final double[][] vertices = deg2rad(new double[][]{
      {0.8537903597029015,-49.60560772150546},
      {359.7264963624239,-49.60560826882517},
      {359.73477207820594,-48.875097433186816},
      {0.8455129889279733,-48.875096899797896}
    });*/
    // Test with negative longitudes
    final double[][] vertices = deg2rad(new double[][]{
      {0.8537903597029015, -49.60560772150546},
      {-0.27350363757605756, -49.60560826882517},
      {-0.2652279217940874, -48.875097433186816},
      {0.8455129889279733, -48.875096899797896}
    });
    final long[] expected = new long[]{9052823,9052827,9052828,9052829,9052830,9052831,9052839,
        9052843,9052844,9052845,9052846,9052847,9052848,9052849,9052850,9052851,9052852,9052853,
        9052854,9052855,9052856,9052857,9052858,9052859,9052860,9052861,9052862,9052863,9052866,
        9052872,9052874,9052896,9052898,9052904,9052906,9054208,9054209,9054210,9054211,9054212,
        9054213,9054214,9054215,9054216,9054217,9054218,9054219,9054220,9054221,9054222,9054223,
        9054224,9054225,9054226,9054227,9054228,9054229,9054230,9054231,9054232,9054233,9054234,
        9054235,9054236,9054237,9054238,9054239,9054240,9054241,9054242,9054243,9054244,9054245,
        9054246,9054247,9054248,9054249,9054250,9054251,9054252,9054253,9054254,9054255,9054256,
        9054257,9054258,9054259,9054260,9054261,9054262,9054263,9054264,9054265,9054266,9054267,
        9054268,9054269,9054270,9054271,9054272,9054274,9054280,9054282,9054304,9054306,9054312,
        9054314,9054336,9054337,9054338,9054339,9054340,9054341,9054342,9054343,9054344,9054345,
        9054346,9054347,9054348,9054349,9054350,9054351,9054352,9054353,9054354,9054355,9054356,
        9054357,9054358,9054359,9054360,9054361,9054362,9054363,9054364,9054368,9054369,9054370,
        9054371,9054372,9054373,9054374,9054375,9054376,9054377,9054378,9054379,9054380,9054384,
        9054400,9054720,11866455,11866461,11867136,11867137,11867138,11867139,11867140,11867141,
        11867142,11867143,11867144,11867145,11867148,11867149,11867152,11867153,11867154,11867155,
        11867156,11867157,11867158,11867159,11867160,11867161,11867164,11867165,11867200,11867201,
        11867202,11867203,11867204,11867205,11867206,11867207,11867208,11867209,11867212,11867213,
        11867216,11867217,11867218,11867219,11867220,11867221,11867222,11867223,11867224,11867225,
        11867228,11867392};
    // - generic algo
    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    // toAladinDraw(depth, bmoc);
    for (int i = 0; it.hasNext(); i++) {
      // System.out.println(it.next());
      assertEquals(expected[i], it.next());
    }
  }
  
  @Test
  public void polyonTest10() {
    final int depth = 12;
    final double[][] vertices = new double[][]{
      { -1.5706045044233712, -0.7295105218483977 },
      { -1.5705168372776197, -0.7295105199399403 },
      { -1.5704291701320918, -0.7295105142145686 },
      { -1.5703415029870114, -0.7295105046722821 },
    };
    final long[] expected = new long[]{195734186, 195734187};
    // - generic algo
    HealpixNested hn = Healpix.getNested(depth);
    HealpixNestedPolygonComputer pc = hn.newPolygonComputer();
    HealpixNestedBMOC bmoc = pc.overlappingCells(vertices);
    assertEquals(expected.length, bmoc.computeDeepSize());
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      assertEquals(expected[i], it.next());
    }
  }
  
  private void toAladinDraw (final int depth, HealpixNestedBMOC bmoc) {
    // NOT STANDARD (SINCE WE DO NOT USE INTERVALS!!)
    System.out.format("draw moc %d/", depth);
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      System.out.format("%d,", it.next());
    }
  }
  
  private static final double[][] deg2rad(final double[][] vertices) {
    for (int i = 0; i < vertices.length; i++) {
      final double[] coos = vertices[i];
      for (int j = 0; j < coos.length; j++) {
        coos[j] = Math.toRadians(coos[j]);
      }
    }
    return vertices;
  }
  
}
