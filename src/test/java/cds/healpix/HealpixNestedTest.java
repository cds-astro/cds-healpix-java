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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.Test;

import cds.healpix.HashComputer;
import cds.healpix.HashParts;
import cds.healpix.Healpix;
import cds.healpix.HealpixNested;
import cds.healpix.HealpixNestedFast;
import cds.healpix.fillingcurve.FillingCurve2DType;

import static cds.healpix.common.math.Math.PI;

import static cds.healpix.Healpix.TRANSITION_LATITUDE;

public final class HealpixNestedTest {

  /*@Test
  public void ijBasePixelDeltaTest() {
    assertEquals(-1, HealpixNestedHashComputer.ijBasePixelDelta( 1,  0)); // d0h = 0
    assertEquals(-1, HealpixNestedHashComputer.ijBasePixelDelta( 2, -1)); // d0h = 1
    assertEquals(-1, HealpixNestedHashComputer.ijBasePixelDelta(-1,  2)); // d0h = 2
    assertEquals(-1, HealpixNestedHashComputer.ijBasePixelDelta( 0,  1)); // d0h = 3
    assertEquals( 0, HealpixNestedHashComputer.ijBasePixelDelta( 0,  0)); // d0h = 4
    assertEquals( 0, HealpixNestedHashComputer.ijBasePixelDelta( 1, -1)); // d0h = 5
    assertEquals( 0, HealpixNestedHashComputer.ijBasePixelDelta( 2, -2)); // d0h = 6
    assertEquals( 0, HealpixNestedHashComputer.ijBasePixelDelta(-2,  2)); // d0h = 6
    assertEquals( 0, HealpixNestedHashComputer.ijBasePixelDelta(-1,  1)); // d0h = 7
    assertEquals( 1, HealpixNestedHashComputer.ijBasePixelDelta( 0, -1)); // d0h = 8
    assertEquals( 1, HealpixNestedHashComputer.ijBasePixelDelta( 1, -2)); // d0h = 9
    assertEquals( 1, HealpixNestedHashComputer.ijBasePixelDelta(-2,  1)); // d0h = 10
    assertEquals( 1, HealpixNestedHashComputer.ijBasePixelDelta(-1,  0)); // d0h = 11
  }*/
  
  @Test
  public void hashTests() throws Exception {
    final int depth = 8;
    HealpixNested hn = Healpix.getNested(depth);
    final long nside = hn.nside;
    final long nsidesq = nside * nside;
    assertEquals(hn.hash(0.03, 0.04), 311384);
    // System.out.println(hn.hash(0.03, 0.04));
    
    assertEquals(nsidesq - 1, hn.hash(Math.toRadians(  0), Math.toRadians(90)));  // 0
    assertEquals(nsidesq - 1, hn.hash(Math.toRadians( 35), Math.toRadians(90)));
    assertEquals(nsidesq - 1, hn.hash(Math.toRadians( 45), Math.toRadians(90))); // pi/4
    assertEquals(nsidesq - 1, hn.hash(Math.toRadians( 55), Math.toRadians(90)));
    assertEquals(2 * nsidesq - 1, hn.hash(Math.toRadians( 90), Math.toRadians(90))); // pi/2
    assertEquals(2 * nsidesq - 1, hn.hash(Math.toRadians(125), Math.toRadians(90)));
    assertEquals(2 * nsidesq - 1, hn.hash(Math.toRadians(135), Math.toRadians(90))); // 3pi/4
    assertEquals(2 * nsidesq - 1, hn.hash(Math.toRadians(155), Math.toRadians(90)));
    assertEquals(3 * nsidesq - 1, hn.hash(Math.toRadians(180), Math.toRadians(90))); // pi
    assertEquals(3 * nsidesq - 1, hn.hash(Math.toRadians(225), Math.toRadians(90))); 
    assertEquals(4 * nsidesq - 1, hn.hash(Math.toRadians(270), Math.toRadians(90))); // 3pi/2 
    assertEquals(4 * nsidesq - 1, hn.hash(Math.toRadians(315), Math.toRadians(90))); 
    assertEquals(nsidesq - 1, hn.hash(Math.toRadians(360), Math.toRadians(90))); // 2pi
    // Test South pole
    assertEquals( 8 * nsidesq, hn.hash(Math.toRadians(  0), Math.toRadians(-90)));  // 0
    assertEquals( 8 * nsidesq, hn.hash(Math.toRadians( 35), Math.toRadians(-90)));
    assertEquals( 8 * nsidesq, hn.hash(Math.toRadians( 45), Math.toRadians(-90))); // pi/4
    assertEquals( 8 * nsidesq, hn.hash(Math.toRadians( 55), Math.toRadians(-90)));
    assertEquals( 9 * nsidesq, hn.hash(Math.toRadians( 90), Math.toRadians(-90))); // pi/2
    assertEquals( 9 * nsidesq, hn.hash(Math.toRadians(125), Math.toRadians(-90)));
    assertEquals( 9 * nsidesq, hn.hash(Math.toRadians(135), Math.toRadians(-90))); // 3pi/4
    assertEquals( 9 * nsidesq, hn.hash(Math.toRadians(155), Math.toRadians(-90)));
    assertEquals(10 * nsidesq, hn.hash(Math.toRadians(180), Math.toRadians(-90))); // pi
    assertEquals(10 * nsidesq, hn.hash(Math.toRadians(225), Math.toRadians(-90))); 
    assertEquals(11 * nsidesq, hn.hash(Math.toRadians(270), Math.toRadians(-90))); // 3pi/2 
    assertEquals(11 * nsidesq, hn.hash(Math.toRadians(315), Math.toRadians(-90))); 
    assertEquals( 8 * nsidesq, hn.hash(Math.toRadians(360), Math.toRadians(-90))); // 2pi
   
    
    HealpixNestedFast hf = new HealpixNestedFast(depth, FillingCurve2DType.Z_ORDER_LUPT);
    
    assertEquals(nsidesq - 1, hf.hash(Math.toRadians(  0), Math.toRadians(90)));  // 0
    assertEquals(nsidesq - 1, hf.hash(Math.toRadians( 35), Math.toRadians(90)));
    assertEquals(nsidesq - 1, hf.hash(Math.toRadians( 45), Math.toRadians(90))); // pi/4
    assertEquals(nsidesq - 1, hf.hash(Math.toRadians( 55), Math.toRadians(90)));
    assertEquals(2 * nsidesq - 1, hf.hash(Math.toRadians( 90), Math.toRadians(90))); // pi/2
    assertEquals(2 * nsidesq - 1, hf.hash(Math.toRadians(125), Math.toRadians(90)));
    assertEquals(2 * nsidesq - 1, hf.hash(Math.toRadians(135), Math.toRadians(90))); // 3pi/4
    assertEquals(2 * nsidesq - 1, hf.hash(Math.toRadians(155), Math.toRadians(90)));
    assertEquals(3 * nsidesq - 1, hf.hash(Math.toRadians(180), Math.toRadians(90))); // pi
    assertEquals(3 * nsidesq - 1, hf.hash(Math.toRadians(225), Math.toRadians(90))); 
    assertEquals(4 * nsidesq - 1, hf.hash(Math.toRadians(270), Math.toRadians(90))); // 3pi/2 
    assertEquals(4 * nsidesq - 1, hf.hash(Math.toRadians(315), Math.toRadians(90))); 
    assertEquals(nsidesq - 1, hf.hash(Math.toRadians(360), Math.toRadians(90))); // 2pi
    // Test South pole
    assertEquals( 8 * nsidesq, hf.hash(Math.toRadians(  0), Math.toRadians(-90)));  // 0
    assertEquals( 8 * nsidesq, hf.hash(Math.toRadians( 35), Math.toRadians(-90)));
    assertEquals( 8 * nsidesq, hf.hash(Math.toRadians( 45), Math.toRadians(-90))); // pi/4
    assertEquals( 8 * nsidesq, hf.hash(Math.toRadians( 55), Math.toRadians(-90)));
    assertEquals( 9 * nsidesq, hf.hash(Math.toRadians( 90), Math.toRadians(-90))); // pi/2
    assertEquals( 9 * nsidesq, hf.hash(Math.toRadians(125), Math.toRadians(-90)));
    assertEquals( 9 * nsidesq, hf.hash(Math.toRadians(135), Math.toRadians(-90))); // 3pi/4
    assertEquals( 9 * nsidesq, hf.hash(Math.toRadians(155), Math.toRadians(-90)));
    assertEquals(10 * nsidesq, hf.hash(Math.toRadians(180), Math.toRadians(-90))); // pi
    assertEquals(10 * nsidesq, hf.hash(Math.toRadians(225), Math.toRadians(-90))); 
    assertEquals(11 * nsidesq, hf.hash(Math.toRadians(270), Math.toRadians(-90))); // 3pi/2 
    assertEquals(11 * nsidesq, hf.hash(Math.toRadians(315), Math.toRadians(-90))); 
    assertEquals( 8 * nsidesq, hf.hash(Math.toRadians(360), Math.toRadians(-90))); // 2pi
    
    
    // System.out.println("Coucou :)");
    
    final HealpixNested hn2 = Healpix.getNested(16);
    HealpixNestedFast hf2 = new HealpixNestedFast(16, FillingCurve2DType.Z_ORDER_LUPT);
    
    // assertEquals(hn2.xMask, hb2.ang2pix(new Pointing(0.5 * PI - 0.7297280099070644, 1.5707963267948966)));
    // assertEquals(hn2.xMask, hn2.hash(1.5707963267948966, 0.7297280099070644));
    // System.out.format(Locale.US, "%.16f +%.16f", Math.toDegrees(1.5707963267948966), Math.toDegrees(0.7297280099070644));
    
    double lb = Math.asin(2.0 / 3.0);
    long h1, h2, h3;
    for (int i = 0; i < 1000000; i++) {
      double lat = lb + (i * 1e-15);
      
      double lon = 2 * PI; //3 * PI / 2;
      // System.out.println("lon: " + Math.toDegrees(lon) + "; lat: " + Math.toDegrees(lat));
      h1 = hn2.hash(lon, lat);
      h2 = hf2.hash(lon, lat);
      
      lon = 0 + 1e-15;
      hn2.hash(lon, lat);
      h1 = hn2.hash(lon, lat);
      h2 = hf2.hash(lon, lat);
      
      lon = 0 - 1e-15;
      hn2.hash(lon, lat);
      h1 = hn2.hash(lon, lat);
      h2 = hf2.hash(lon, lat);
    }
  
    
    /*double lat = +PI * 0.4852;
    long signLat = toBits(lat);                             
    double absLat = fromBits(signLat & BUT_SIGN_BIT_MASK_L); 
    signLat &= SIGN_BIT_MASK_L;
    */
    /*int zeroOr8 = (int) (signLat >> 60);
    System.out.println((1 >> zeroOr8));
    System.out.println("----");
    System.out.println("-0: " + ((int) -0));
    System.out.println("-1: " + ((int) -1 ));
    System.out.println("-2: " + ((int) -2));
    System.out.println("-1.25: " + ((int) -1.25));
    System.out.println("-2.75: " + ((int) -2.75));
    System.out.println("----");
    System.out.println(signLat);
    System.out.println((signLat >>> 60));
    System.out.println((2 >>> (signLat >>> 60)));*/
    
  }

  
  @Test 
  public void noNeighbourValueTest() {
    HealpixNested hn = new HealpixNested(5, FillingCurve2DType.Z_ORDER_LUPT);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    hn = new HealpixNested(5, FillingCurve2DType.Z_ORDER_XOR);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    hn = new HealpixNested(5, FillingCurve2DType.Z_ORDER_OR);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    //
    hn = new HealpixNested(15, FillingCurve2DType.Z_ORDER_LUPT);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    hn = new HealpixNested(15, FillingCurve2DType.Z_ORDER_XOR);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    hn = new HealpixNested(15, FillingCurve2DType.Z_ORDER_OR);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    //
    hn = new HealpixNested(25, FillingCurve2DType.Z_ORDER_LUPT);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    hn = new HealpixNested(25, FillingCurve2DType.Z_ORDER_XOR);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
    hn = new HealpixNested(25, FillingCurve2DType.Z_ORDER_OR);
    assertEquals(-1, hn.encodeRegularHash(-1, -1, -1));
  }
  
  @Test
  public void regularRingEncodeTest() throws IOException {
    final int depth = 7;
    final HealpixNested hn = Healpix.getNested(depth);
    final List<String> lines = new ArrayList<String>((int) hn.nHash + 1);
    lines.add("nested,ring,ra,dec");
    for (long i = 0; i < hn.nHash; i++) {
      // System.out.println("hashNested: " + i);
      final HashParts hparts = hn.decodeRegularHash(i);
      final double[] center = hn.newVerticesAndPathComputer().center(i);  
      long ring = hn.regularRingEncode(hparts.baseCellHash(), hparts.iInBaseCell(), hparts.jInBaseCell());
      // lines.add(i + "," + ring + "," + String.format(Locale.US, "%.6f,%.6f", center[0], center[1]));
      final HashParts decodedRing = hn.decodeRegularRing(ring);
      if (decodedRing.baseCellHash() != hparts.baseCellHash()
          || decodedRing.iInBaseCell() != hparts.iInBaseCell()
          || decodedRing.jInBaseCell() != hparts.jInBaseCell()) {
        throw new Error("Decode error. Expected: " + hparts + "; Actual: " + decodedRing);
      }
    }
    // final File f = new File("./nested2ring.depth" + depth + ".csv");
    // Files.write(f.toPath(), lines, StandardCharsets.US_ASCII);
  }

  
}
