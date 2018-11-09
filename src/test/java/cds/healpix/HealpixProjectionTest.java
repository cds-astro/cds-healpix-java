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


import static cds.healpix.common.math.Math.atan2;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.sqrt;
import static cds.healpix.common.math.Math.toDegrees;
import static cds.healpix.common.math.Math.toRadians;

import java.io.File;
import java.io.IOException;
/*import java.nio.charset.StandardCharsets;
import java.nio.file.Files;*/
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.junit.Test;

import cds.healpix.HealpixProjection;

public class HealpixProjectionTest {

  private final HealpixProjection HPX = new HealpixProjection();
  
  
  public static final class Point {

    private double lon, lat; // RA, Dec in radians
    private double[] xy = new double[2];
    private double color;

    public Point(final double lon, final double lat) {
      this(lon, lat, (90 + toDegrees(lat)) / 180.0);
    }

    public Point(final double lon, final double lat, final double color) {
      this.lon = lon;
      this.lat = lat;
      this.color = color;
    }

    public void setLonLat(
        final double longitudeRadians, final double latitudeRadians) {
      this.lon = longitudeRadians;
      this.lat = latitudeRadians;
      this.color = (90 + toDegrees(lat)) / 180.0;
    }

    @Override public String toString() {
      return String.format(Locale.US, "%.6f,%.6f,%.6f,%.6f,%.4f", 
          toDegrees(lon), toDegrees(lat), xy[0], xy[1], color);
    }

  }

  public static List<Point> buildPoints() {
    final List<Point> p = new ArrayList<Point>(1800 * 12 + 3600 * 11);
    // Meridians (lon = cte)
    for (int l = 0; l <= 360; l += 30) { // 0, 30, 60, 90, ... 360 // (+30 x12)
      for (int b = -900; b <= +900; b++) { // 180 * 10 = 1800 points
        p.add(new Point(toRadians(l), toRadians(b / 10.0)));
        if (l > 0 && l < 360) {
          p.add(new Point(toRadians(l) + 1.e-14, toRadians(b / 10.0)));
          p.add(new Point(toRadians(l) - 1.e-14, toRadians(b / 10.0)));
        } else if (l == 360) {
          p.add(new Point(toRadians(l) - 1.e-14, toRadians(b / 10.0)));
        }
      }
    }
    // Circles of latitude (lat = cte)
    for (int b = -90; b <= +90; b += 15) { // +- 0, 15, 30, 45, 60, 75 // (+-15 1+2x5) 
      for (int l = 0; l < 3600; l++) { 
        p.add(new Point(toRadians(l / 10.0), toRadians(b)));
      }
    }
    // lat = 90 - epislon
    for (int l = 0; l < 3600; l++) { 
      p.add(new Point(toRadians(l / 10.0), toRadians(90) - 1.e-12));
    }
    // lat = -90 + epsilon
    for (int l = 0; l < 3600; l++) { 
      p.add(new Point(toRadians(l / 10.0), toRadians(-90) + 1.e-12));
    }
    // Add center of HEALPix cells of level 3
    /*final HealpixNested hn = new HealpixNested(3);
    for (long h = 0; h < hn.nHashes(); h++) {
        final cds.common.astrometry.LonLat ll = hn.center(h);
        p.add(new Point(ll.lon(), ll.lat(), -1));
    }
    // Add vertices of HEALPix cells of level 3
    for (long h = 0; h < hn.nHashes(); h++) {
        for (final cds.common.astrometry.LonLat ll: hn.vertices(h)) {
            p.add(new Point(ll.lon(), ll.lat(), -2));
        }
    }*/
    return p;
  }

  public static final double vincentyDist(
      final double lon1Rad, final double lat1Rad,
      final double lon2Rad, final double lat2Rad) {
  final double deltaLon = lon2Rad - lon1Rad;
  final double cDeltaLon = cos(deltaLon);
  final double clat1 = cos(lat1Rad); 
  final double slat1 = sin(lat1Rad);
  final double clat2 = cos(lat2Rad); 
  final double slat2 = sin(lat2Rad);
  final double a = clat2 * sin(deltaLon);
  final double b = clat1 * slat2 - slat1 * clat2 * cDeltaLon;
  return atan2(sqrt(a * a + b * b),
          slat1 * slat2 + clat1 * clat2 * cDeltaLon);
}
  
  @Test
  public void testProj() throws IOException {
    // Build points
    final List<Point> points = buildPoints();
    // Do projections
    final List<String> s = new ArrayList<String>(points.size() + 1);
    s.add("lon,lat,x,y,color");
    final StringBuilder strBuild = new StringBuilder();
    int n = 0;
    for (final Point p : points) {
        strBuild.delete(0, strBuild.length());
        HPX.project(p.lon, p.lat, p.xy);
        double[] lonlat = HPX.unproject(p.xy[0], p.xy[1]);
        double d = vincentyDist(p.lon, p.lat, lonlat[0], lonlat[1]);
        assertFalse(d > toRadians(1 / (3600 * 1000 * 1000d)));
        /*if (d > toRadians(1 / (3600 * 1000 * 1000d))) {
          throw new Error("Oups: " + toDegrees(d * 3600) + "; " + p.lon + "; " + p.lat + "; x: " + p.xy[0] + "; y: " + p.xy[1] + "; " + lonlat[0] + "; " + lonlat[1]);
        }*/
        assertEquals(p.lat, lonlat[1], toRadians(1 / (3600 * 1000 * 1000d)));
        /*if (Math.abs(p.lat - lonlat[1]) > toRadians(1 / (3600 * 1000 * 1000d))) {
          throw new Error("lat: " + p.lat + " != " + lonlat[1]);
        }*/
        
        
        /*if (Math.abs(p.lat) != 0.5 * Math.PI && Math.abs(p.lon - lonlat[0]) > toRadians(1)) {
          System.out.println("l: " + p.lon + "; b:  " + p.lat + "; x:" + p.xy[0] + "; y: " + p.xy[1] + "; rl: " + lonlat[0] + "; rlat: " + lonlat[1]);
          n++;
        }*/
        // System.out.println("ok");
        // HealpixUtil.proj(ll, point);
        strBuild.append(p.toString());
        // To string to build output file
        s.add(strBuild.toString());
    }
    // System.out.println("N = " + n);
    // Write result file
    // Files.write(new File("test.healpix.proj.csv").toPath(), s, StandardCharsets.US_ASCII);
  }
  
  /*public static void main(final String[] args) throws IOException {
    final HealpixProjectionTest t = new HealpixProjectionTest();
    t.testProj();
  }*/

  
}
