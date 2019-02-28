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

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import cds.healpix.Healpix;
import static cds.healpix.Healpix.*;
import static cds.healpix.common.math.Math.atan2;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.sqrt;
import static cds.healpix.common.math.Math.toDegrees;
import static cds.healpix.common.math.Math.toRadians;
import static cds.healpix.VerticesAndPathComputer.ALL_CARDINAL_POINTS;
import static org.junit.Assert.*;

import java.util.List;
import java.util.ArrayList;
/*import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Locale;*/

public final class HealpixTest {

  @Rule
  public ExpectedException exception = ExpectedException.none();

  @Test
  public void testConstants() {
    assertEquals(DEPTH_MAX, 29);
    assertEquals(TRANSITION_Z, 2 / 3d, 0);
    assertEquals(TRANSITION_LATITUDE, Math.asin(TRANSITION_Z), 0);
  }

  @Test
  public void nsideTest() {
    assertEquals(1, nside(0));
    assertEquals(2, nside(1));
    assertEquals(4, nside(2));
    assertEquals(8, nside(3));
    assertEquals(16, nside(4));
    assertEquals(32, nside(5));
    assertEquals(64, nside(6));
    assertEquals(128, nside(7));
    assertEquals(256, nside(8));
    assertEquals(512, nside(9));
    assertEquals(1024, nside(10));
    assertEquals(2048, nside(11));
    assertEquals(4096, nside(12));
    assertEquals(8192, nside(13));
    assertEquals(16384, nside(14));
    assertEquals(32768, nside(15));
    assertEquals(65536, nside(16));
    assertEquals(131072, nside(17));
    assertEquals(262144, nside(18));
    assertEquals(524288, nside(19));
    assertEquals(1048576, nside(20));
    assertEquals(2097152, nside(21));
    assertEquals(4194304, nside(22));
    assertEquals(8388608, nside(23));
    assertEquals(16777216, nside(24));
    assertEquals(33554432, nside(25));
    assertEquals(67108864, nside(26));
    assertEquals(134217728, nside(27));
    assertEquals(268435456, nside(28));
    assertEquals(536870912, nside(29));
  }
  @Test
  public void nsideTestExceptionMin1() {
    exception.expect(IllegalArgumentException.class);
    exception.expectMessage("Expected depth in [0, 29]. Actual: -1");
    nside(-1);
  }
  @Test
  public void nsideTestException30() {
    exception.expect(IllegalArgumentException.class);
    exception.expectMessage("Expected depth in [0, 29]. Actual: 30");
    nside(30);
  }

  @Test
  public void depthTest() {
    assertEquals(0, depth(1));
    assertEquals(1, depth(2));
    assertEquals(2, depth(4));
    assertEquals(3, depth(8));
    assertEquals(4, depth(16));
    assertEquals(5, depth(32));
    assertEquals(6, depth(64));
    assertEquals(7, depth(128));
    assertEquals(8, depth(256));
    assertEquals(9, depth(512));
    assertEquals(10, depth(1024));
    assertEquals(11, depth(2048));
    assertEquals(12, depth(4096));
    assertEquals(13, depth(8192));
    assertEquals(14, depth(16384));
    assertEquals(15, depth(32768));
    assertEquals(16, depth(65536));
    assertEquals(17, depth(131072));
    assertEquals(18, depth(262144));
    assertEquals(19, depth(524288));
    assertEquals(20, depth(1048576));
    assertEquals(21, depth(2097152));
    assertEquals(22, depth(4194304));
    assertEquals(23, depth(8388608));
    assertEquals(24, depth(16777216));
    assertEquals(25, depth(33554432));
    assertEquals(26, depth(67108864));
    assertEquals(27, depth(134217728));
    assertEquals(28, depth(268435456));
    assertEquals(29, depth(536870912));
  }

  @Test
  public void depthTest0() {
    exception.expect(IllegalArgumentException.class);
    exception.expectMessage("Nside must be a power of 2 in [1-2^29]");
    depth(0);
  }
  @Test
  public void depthTest3() {
    exception.expect(IllegalArgumentException.class);
    exception.expectMessage("Nside must be a power of 2 in [1-2^29]");
    depth(3);
  }

  @Test
  public void crossTestNsideDepth() {
    for (int depth = 0; depth <= Healpix.DEPTH_MAX; depth++) {
      assertEquals(depth, depth(nside(depth)));
    }
  }

  @Test
  public void nHashTest() {
    assertEquals(12L, nHash(0));
    assertEquals(48L, nHash(1));
    assertEquals(192L, nHash(2));
    assertEquals(768L, nHash(3));
    assertEquals(3072L, nHash(4));
    assertEquals(12288L, nHash(5));
    assertEquals(49152L, nHash(6));
    assertEquals(196608L, nHash(7));
    assertEquals(786432L, nHash(8));
    assertEquals(3145728L, nHash(9));
    assertEquals(12582912L, nHash(10));
    assertEquals(50331648L, nHash(11));
    assertEquals(201326592L, nHash(12));
    assertEquals(805306368L, nHash(13));
    assertEquals(3221225472L, nHash(14));
    assertEquals(12884901888L, nHash(15));
    assertEquals(51539607552L, nHash(16));
    assertEquals(206158430208L, nHash(17));
    assertEquals(824633720832L, nHash(18));
    assertEquals(3298534883328L, nHash(19));
    assertEquals(13194139533312L, nHash(20));
    assertEquals(52776558133248L, nHash(21));
    assertEquals(211106232532992L, nHash(22));
    assertEquals(844424930131968L, nHash(23));
    assertEquals(3377699720527872L, nHash(24));
    assertEquals(13510798882111488L, nHash(25));
    assertEquals(54043195528445952L, nHash(26));
    assertEquals(216172782113783808L, nHash(27));
    assertEquals(864691128455135232L, nHash(28));
    assertEquals(3458764513820540928L, nHash(29));
  }

  @Test
  public void nIsolatitudeTest() {
    assertEquals(3, Healpix.nIsolatitudeRings(0));
    assertEquals(7, Healpix.nIsolatitudeRings(1));
    assertEquals(15, Healpix.nIsolatitudeRings(2));
    assertEquals(31, Healpix.nIsolatitudeRings(3));
    assertEquals(63, Healpix.nIsolatitudeRings(4));
    assertEquals(127, Healpix.nIsolatitudeRings(5));
    assertEquals(255, Healpix.nIsolatitudeRings(6));
    assertEquals(511, Healpix.nIsolatitudeRings(7));
    assertEquals(1023, Healpix.nIsolatitudeRings(8));
    assertEquals(2047, Healpix.nIsolatitudeRings(9));
    assertEquals(4095, Healpix.nIsolatitudeRings(10));
    assertEquals(8191, Healpix.nIsolatitudeRings(11));
    assertEquals(16383, Healpix.nIsolatitudeRings(12));
    assertEquals(32767, Healpix.nIsolatitudeRings(13));
    assertEquals(65535, Healpix.nIsolatitudeRings(14));
    assertEquals(131071, Healpix.nIsolatitudeRings(15));
    assertEquals(262143, Healpix.nIsolatitudeRings(16));
    assertEquals(524287, Healpix.nIsolatitudeRings(17));
    assertEquals(1048575, Healpix.nIsolatitudeRings(18));
    assertEquals(2097151, Healpix.nIsolatitudeRings(19));
    assertEquals(4194303, Healpix.nIsolatitudeRings(20));
    assertEquals(8388607, Healpix.nIsolatitudeRings(21));
    assertEquals(16777215, Healpix.nIsolatitudeRings(22));
    assertEquals(33554431, Healpix.nIsolatitudeRings(23));
    assertEquals(67108863, Healpix.nIsolatitudeRings(24));
    assertEquals(134217727, Healpix.nIsolatitudeRings(25));
    assertEquals(268435455, Healpix.nIsolatitudeRings(26));
    assertEquals(536870911, Healpix.nIsolatitudeRings(27));
    assertEquals(1073741823, Healpix.nIsolatitudeRings(28));
    assertEquals(2147483647, Healpix.nIsolatitudeRings(29));
  }

  @Test
  public void nHashTestMin1() {
    exception.expect(IllegalArgumentException.class);
    exception.expectMessage("Expected depth in [0, 29]. Actual: -1");
    nHash(-1);
  }

  @Test
  public void nHashTestMin30() {
    exception.expect(IllegalArgumentException.class);
    exception.expectMessage("Expected depth in [0, 29]. Actual: 30");
    nHash(30);
  }

  private static List<double[]> buildPoints() {
    final List<double[]> p = new ArrayList<double[]>(1800 * 12 + 3600 * 11);
    // Meridians (lon = cte)
    for (int l = 0; l <= 360; l += 30) { // 0, 30, 60, 90, ... 360 // (+30 x12)
      for (int b = -900; b <= +900; b++) { // 180 * 10 = 1800 points
        p.add(new double[]{toRadians(l), toRadians(b / 10.0)});
        if (l > 0 && l < 360) {
          p.add(new double[]{toRadians(l) + 1.e-14, toRadians(b / 10.0)});
          p.add(new double[]{toRadians(l) - 1.e-14, toRadians(b / 10.0)});
        } else if (l == 360) {
          p.add(new double[]{toRadians(l) - 1.e-14, toRadians(b / 10.0)});
        }
      }
    }
    // Circles of latitude (lat = cte)
    for (int b = -90; b <= +90; b += 15) { // +- 0, 15, 30, 45, 60, 75 // (+-15 1+2x5) 
      for (int l = 0; l < 3600; l++) { 
        p.add(new double[]{toRadians(l / 10.0), toRadians(b)});
      }
    }
    // lat = 90 - epislon
    for (int l = 0; l < 3600; l++) { 
      p.add(new double[]{toRadians(l / 10.0), toRadians(90) - 1.e-12});
    }
    // lat = -90 + epsilon
    for (int l = 0; l < 3600; l++) { 
      p.add(new double[]{toRadians(l / 10.0), toRadians(-90) + 1.e-12});
    }
    // Add center and vertices of HEALPix cells of level 3
    final HealpixNested hn = Healpix.getNested(3);
    for (long h = 0; h < hn.nHash; h++) {
      p.add(hn.center(h));
      for (double[] ll : hn.vertices(h, ALL_CARDINAL_POINTS).values()) {
        p.add(ll);
      }
    }
    return p;
  }

  private static final double vincentyDist(
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

  /* TEST REMOVED TO BE COMPATIBLE WITH JAVA 6
  @Test
  public void testProj() throws IOException {
    // Build points
    final List<double[]> points = buildPoints();
    // Do projections
    final List<String> s = new ArrayList<String>(points.size() + 1);
    s.add("orgLonDeg,orgLatDeg,projX,projY,unpojLonDeg,unprojLatDeg");
    for (final double[] point : points) {
      final double[] xy = Healpix.UI.project(point[LON_INDEX], point[LAT_INDEX]);
      final double[] radec = Healpix.UI.unproject(xy[X_INDEX], xy[Y_INDEX]);
      double d = vincentyDist(point[LON_INDEX], point[LAT_INDEX], radec[LON_INDEX], radec[LAT_INDEX]);
      if (d > toRadians(1 / (3600 * 1000 * 1000d))) { // d > 1 microarcsec
        throw new Error("\nOrg. lon: " + toDegrees(point[LON_INDEX]) + "; lat: " + toDegrees(point[LAT_INDEX])
            + "\nProj. x: " + xy[X_INDEX] + "; y: " + xy[Y_INDEX]
            + "\nDeproj. lon: " + toDegrees(radec[LON_INDEX]) + "; lat: " + toDegrees(radec[LAT_INDEX])
            + "\nAngularDist: " + toDegrees(d * 3600) + " arcsec. " + s.size());
      }
      // To string to build output file
      s.add(String.format(Locale.US, "%.13f,%.13f,%.13f,%.13f,%.13f,%.13f",
          toDegrees(point[LON_INDEX]), toDegrees(point[LAT_INDEX]), xy[X_INDEX], xy[Y_INDEX],
          toDegrees(radec[LON_INDEX]), toDegrees(radec[LAT_INDEX])));
    }
    // Write result file
    Files.write(new File("target/test-results/healpix.proj.actual.csv").toPath(), s, StandardCharsets.US_ASCII);
    byte[] f1 = Files.readAllBytes(new File("target/test-results/healpix.proj.actual.csv").toPath());
    byte[] f2 = Files.readAllBytes(new File("src/test/resources/healpix.proj.expected.csv").toPath());
    assertTrue(Arrays.equals(f1,  f2));
  }*/

  @Test
  public void uniqTest() {
    int depth = 10;
    long nHash = Healpix.nHash(depth);
    for (long h = 0; h < nHash; h++) {
      long u = Healpix.uniq(depth, h);
      assertEquals(depth, Healpix.uniq2depth(u));
      assertEquals(h, Healpix.uniq2hash(u));
      assertEquals(h, Healpix.uniq2hash(u, depth));
    }
  }
  
}
