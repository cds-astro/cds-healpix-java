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

  @Test
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
        //"287.18290328,-02.26765169",
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
  }
  
}
