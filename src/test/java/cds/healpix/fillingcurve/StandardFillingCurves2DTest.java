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

import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_BYTE;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_INT;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_SHORT;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_OR_BYTE;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_OR_INT;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_OR_SHORT;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_XOR_BYTE;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_XOR_INT;
import static cds.healpix.fillingcurve.ZOrderCurve2DImpls.ZOC_VMSB_XOR_SHORT;
import static org.junit.Assert.*;

import java.util.Random;

import org.junit.Test;


public final class StandardFillingCurves2DTest {

  @Test
  public void testIJ2HashByte() {
    for (int i = 0; i < 256; i++) {
      for (int j = 0; j < 256; j++) {
        long h1 = ZOC_VMSB_OR_BYTE.ij2hash(i, j);
        long h2 = ZOC_VMSB_XOR_BYTE.ij2hash(i, j);
        long h3 = ZOC_VMSB_LOOKUP_BYTE.ij2hash(i, j);
        assertEquals(h1, h2);
        assertEquals(h1, h3);
      }
    }
  }

  @Test
  public void testIJ2HashShort() {
    final long max = (1 << 16) - 1;
    for (int i = 0; i <= max; i++) {
      for (int j = 0; j <= max; j++) {
        long h1 = ZOC_VMSB_OR_SHORT.ij2hash(i, j);
        long h2 = ZOC_VMSB_XOR_SHORT.ij2hash(i, j);
        long h3 = ZOC_VMSB_LOOKUP_SHORT.ij2hash(i, j);
        assertEquals(h1, h2);
        assertEquals(h1, h3);
      }
    }
  }
  
  @Test
  public void testIJ2HashInt() {
    final long max = 20000;
    final Random r = new Random();
    for (int i = 0; i <= max; i++) {
      int k = Math.abs(r.nextInt());
      for (int j = 0; j <= max; j++) {
        int l = Math.abs(r.nextInt());
        long h1 = ZOC_VMSB_OR_INT.ij2hash(k, l);
        long h2 = ZOC_VMSB_XOR_INT.ij2hash(k, l);
        long h3 = ZOC_VMSB_LOOKUP_INT.ij2hash(k, l);
        assertEquals(h1, h2);
        assertEquals(h1, h3);
      }
    }
  }
  
  @Test
  public void testHash2ijByte() {
    final long max = (1 << 16) - 1;
    for (long h = 0; h <= max; h++) {
        long ij1 = ZOC_VMSB_OR_BYTE.hash2ij(h);
        long ij2 = ZOC_VMSB_XOR_BYTE.hash2ij(h);
        long ij3 = ZOC_VMSB_LOOKUP_BYTE.hash2ij(h);
        assertEquals(ij1, ij2);
        assertEquals(ij1, ij3);
    }
  }
  
  @Test
  public void testHash2ijShort() {
    final long max = (1L << 32) - 1;
    System.out.println("max: " + max);
    for (long h = 0; h <= max; h++) {
        long ij1 = ZOC_VMSB_OR_SHORT.hash2ij(h);
        long ij2 = ZOC_VMSB_XOR_SHORT.hash2ij(h);
        long ij3 = ZOC_VMSB_LOOKUP_SHORT.hash2ij(h);
        assertEquals(ij1, ij2);
        assertEquals(ij1, ij3);
    }
  }
  
  @Test
  public void testHash2ijInt() {
    final long max = 20000000;
    System.out.println("max: " + max);
    final Random r = new Random();
    for (long h = 0; h <= max; h++) {
        final long rdm = Math.abs(r.nextLong());
        long ij1 = ZOC_VMSB_OR_INT.hash2ij(rdm);
        long ij2 = ZOC_VMSB_XOR_INT.hash2ij(rdm);
        long ij3 = ZOC_VMSB_LOOKUP_INT.hash2ij(rdm);
        assertEquals(ij1, ij2);
        assertEquals(ij1, ij3);
    }
  }
  
}
