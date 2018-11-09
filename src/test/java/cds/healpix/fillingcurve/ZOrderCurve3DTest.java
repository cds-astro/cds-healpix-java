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

package cds.healpix.fillingcurve;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import cds.healpix.fillingcurve.ZOrderCurve3D;

public class ZOrderCurve3DTest {
  
  @Test
  public void testConv() {
    for (int n = 0; n < 10000; n++) {
      int i = (int) (Math.random() * 524288);
      int j = (int) (Math.random() * 524288);
      int k = (int) (Math.random() * 524288);
      long hash = ZOrderCurve3D.INSTANCE.ijk2hash(i, j, k);
      //System.out.println("i:" + i + "; j: " + j + "; k: " + k);
      long ijk1 = (((long) k) << 40) | (((long) j) << 20) | i;
      /*System.out.println("ijk:" + ijk1);
      System.out.println("hash:" + hash);*/
      long ijk = ZOrderCurve3D.INSTANCE.hash2ijk(hash);
      /*System.out.println("ijk:" + ijk);
      System.out.format("%19s\n", Integer.toBinaryString(i));
      System.out.format("%19s\n", Integer.toBinaryString(j));
      System.out.format("%19s\n", Integer.toBinaryString(k));
      System.out.format("%57s\n", Long.toBinaryString(hash));*/
      assertEquals(i, ZOrderCurve3D.INSTANCE.ijk2i(ijk));
      assertEquals(j, ZOrderCurve3D.INSTANCE.ijk2j(ijk));
      assertEquals(k, ZOrderCurve3D.INSTANCE.ijk2k(ijk));

    }
  }

  @Test
  public void testLaceUnlace() {
    for (int i = 0; i <= 255; i++) {
      final int t = (int) lace(0, 0, i);
      final int ip = (int) unlace((long) t);
      assertEquals(i, ip);
    }
  }

  @Test
  public void testLaceUnlace2() {
    for (int i = 0; i <= 7; i++) {
      for (int j = 0; j <= 7; j++) {
        for (int k = 0; k <= 7; k++) {
          final long ijk = ((long) k) << 40 | ((long) j) << 20 | i;
          final int t = (int) lace(k, j, i);
          final long ip = unlace((long) t);
          assertEquals(ijk, ip);
        }
      }
    }
  }


  private static final long lace(final int i, final int j, final int k) {
    // We consider that i, j and k elements are bytes!!
    int l = (i << 16) | (j << 8) | k;
    l = (l & 0x000000F0) << 8 | (l & 0x000F0000) >> 8 | (l & 0x0000F000) << 4 | (l & 0x00000F00) >> 4 | l & 0x00F0000F;
    l = (l & 0x0000C00C) << 4 | (l & 0x00300300) >> 4 | (l & 0x000C00C0) << 2 | (l & 0x00030030) >> 2 | l & 0x00C03C03;
    l = (l & 0x00082082) << 2 | (l & 0x00410410) >> 2 | (l & 0x00208208) << 1 | (l & 0x00104104) >> 1 | l & 0x00861861;
    return l;
  }

  private static long unlace(long hash) {
    int l = (int) hash;
    l = (l & 0x00082082) << 1 | (l & 0x00410410) >> 1 | (l & 0x00208208) >> 2 | (l & 0x00104104) << 2 | l & 0x00861861;
    l = (l & 0x0000C00C) << 2 | (l & 0x00300300) >> 2 | (l & 0x000C00C0) >> 4 | (l & 0x00030030) << 4 | l & 0x00C03C03;
    l = (l & 0x000000F0) << 4 | (l & 0x000F0000) >> 4 | (l & 0x0000F000) >> 8 | (l & 0x00000F00) << 8 | l & 0x00F0000F;
    return (((long) (l >> 16)) << 40) | (((long) (l >> 8) & 0x000000FF) << 20) | (l & 0x000000FF);
  }

  private static final void buildLUPT_TO_HASH() {
    System.out.println("UPT_TO_HASH values:");
    for (int i = 0; i <= 255; i++) {
      final int t = (int) lace(0, 0, i);
      // System.out.format("0x%05X --> %16s\n", t, Integer.toBinaryString(t).replace(' ', '0'));
      System.out.format("0x%06X, ", t); // Integer.toHexString(t)
    }
  }
  
  private static final void buildLUPT_TO_IJK_INT() {
    System.out.println("LUPT_TO_IJK_INT values:");
    for (int i = 0; i < 512; i++) {
      final long t = unlace(i);
      // System.out.format(" %8s: 0x%011X --> %64s\n", Long.toBinaryString(i).replace(' ', '0'), t, Long.toBinaryString(t).replace(' ', '0'));
      System.out.format("0x%011XL, ", t); // Integer.toHexString(t)
    }
  }
  
  public static void main(final String[] args) {
    buildLUPT_TO_HASH();
    buildLUPT_TO_IJK_INT();
  }
  
}
