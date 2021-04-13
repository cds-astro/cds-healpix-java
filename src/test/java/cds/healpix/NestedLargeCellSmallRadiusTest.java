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

import java.io.IOException;

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.HealpixNestedBMOC;
import cds.healpix.NestedLargeCellSmallRadius;

public class NestedLargeCellSmallRadiusTest {

  // First test coherency with NestedLageCellTest

  public void coneTest(double coneCenterLonDeg, double coneCenterLatDeg, double radiusDeg, int depth,
      final long[] expectedRes) {
    System.err.println("Cone test. Center "
        + "lon: "  + coneCenterLonDeg + "deg, "
        + "lat: "  + coneCenterLatDeg + "deg, "
        + "radius: " + radiusDeg + " deg, "
        + "depth: " + depth
        );
    final double coneCenterLonRad = Math.toRadians(coneCenterLonDeg);
    final double coneCenterLatRad = Math.toRadians(coneCenterLatDeg);
    final double radiusRad = Math.toRadians(radiusDeg);
    final NestedLargeCellSmallRadius nlc = new NestedLargeCellSmallRadius(Healpix.getNested(depth), radiusRad);
    final HealpixNestedBMOC moc = nlc.overlappingCells(coneCenterLonRad, coneCenterLatRad);
    int i = 0;
    for (final HealpixNestedBMOC.CurrentValueAccessor cell : moc) {
      // System.out.println(cell);
      assertEquals(expectedRes[i++], cell.getHash());
    }
    assertEquals(expectedRes.length, i);
  }

  @Test
  public void coneTestEquatorialZone() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw cirlcle(076.30808, +08.25603, 7.601')
    coneTest(076.30808, +08.25603, 7.601 / 60d, 8, new long[]{374034, 374040, 374041, 374042});
    // LogManager manager = LogManager.getLogManager();
    // manager.readConfiguration();
  }

  @Test
  public void coneTestNorthPolarCap1() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(117.82775, +75.34039, 1.996 deg)
    coneTest(117.82775, +75.34039, 1.996, 4, new long[]{493, 498, 504, 505});
  }

  @Test
  public void coneTestNorthPolarCap2() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(119.27917, +75.66886, 2.692 deg)
    coneTest(119.27917, +75.66886, 2.692, 4, new long[]{493, 498, 499, 504, 505, 506});

  }

  @Test
  public void coneTestNorthPolarCap3() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(055.75724, +89.20594, 37.04')
    coneTest(055.75724, +89.20594, 37.04 / 60d, 6, new long[]{4093, 4094, 4095, 8190, 8191});
  }

  @Test
  public void coneTestNorthPolarCap4() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(089.91417, +47.04352, 59')
    coneTest(089.91417, +47.04352, 59 / 60d, 5, new long[]{349, 351, 373, 1710, 1711, 1722});
  }

  @Test
  public void coneTestNorthPolarCap5() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(086.88067, +45.73794, 40.23')
    coneTest(086.88067, +45.73794, 40.23 / 60d, 5, new long[]{349, 350, 351, 372});
  }

  @Test
  public void coneTestNorthPolarCap6() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(280.75724, +89.20594, 54.04')
    coneTest(280.75724, +89.20594, 54.04 / 60d, 5, new long[]{1023, 2047, 3069, 3071, 4094, 4095});
  }

  @Test
  public void coneTestNorthPolarCap7() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(149.20825, +89.36753, 0.600 arcsec)
    coneTest(149.20825, +89.36753, 0.6 / 3600d, 18, new long[]{137429366799L, 137429366810L
        , 137429366820L, 137429366821L, 137429366822L, 137429366823L, 137429366832L});
  }


  @Test
  public void coneTestNorthPolarCap8() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(134.41234, +89.99997, 0.085'')
    coneTest(134.41234, +89.99997, 0.085 / 3600d, 19,
        new long[]{274877906943L, 549755813887L, 824633720831L});
  }


  @Test
  public void coneTestSouthPolarCap1() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(117.82775, -75.34039, 1.996 deg)
    coneTest(117.82775, -75.34039, 1.996, 4, new long[]{2313, 2315, 2318, 2337});
  }

  @Test
  public void coneTestEquatorialAndSouthPolarCap2() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(89.50838, +42.21108, 53.72')
    coneTest(89.50838, +42.21108, 53.72 / 60d, 5, new long[]{341, 343, 1706, 6143});
  }

  @Test
  public void coneTestEquatorialAndSouthPolarCap3() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(88.57, -41.86, 59')
    coneTest(88.57, -41.86, 59.0 / 60d, 5, new long[]{5120, 8532, 8533, 8535});
  }

  @Test
  public void coneTestEquatorialAndSouthPolarCap4() throws SecurityException, IOException {
    System.setProperty("java.util.logging.config.file", "src/test/resources/logging.properties");
    // draw circle(088.71829, -41.52463, 59')
    coneTest(88.71829, -41.52463, 59.0 / 60d, 5, new long[]{5120, 8532, 8533, 8535});
  }
}
