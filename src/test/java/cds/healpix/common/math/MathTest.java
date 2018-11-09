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

package cds.healpix.common.math;

import org.junit.Test;

import static cds.healpix.common.math.Math.*;
import static org.junit.Assert.*;

public final class MathTest {


  @Test
  public void oneMinusSinTest() {
    oneMinusSinSingleTest(0);
    oneMinusSinSingleTest(HALF_PI);
    oneMinusSinSingleTest(-HALF_PI);
    oneMinusSinSingleTest(PI);
    oneMinusSinSingleTest(-PI);
    oneMinusSinSingleTest(TWO_PI);
    oneMinusSinSingleTest(-TWO_PI);
    for (int i = 0; i < 10000; i++) {
      final double a = java.lang.Math.random() * TWO_PI;
      oneMinusSinSingleTest(a);
      oneMinusSinSingleTest(-a);
    }
  }
  
  private static final void oneMinusSinSingleTest(final double a) {
    // System.out.println("aRad: " + a);
    double res = oneMinusSin(a);
    assertEquals(1 - sin(a), res, 1e-15);
    assertTrue(res >= 0);
    assertTrue(res <= 2);
  }
  
  @Test
  public void sqrtOfOneMinusSinTest() {
    sqrtOneMinusSinSingleTest(0);
    sqrtOneMinusSinSingleTest(HALF_PI);
    sqrtOneMinusSinSingleTest(-HALF_PI);
    sqrtOneMinusSinSingleTest(PI);
    sqrtOneMinusSinSingleTest(-PI);
    sqrtOneMinusSinSingleTest(TWO_PI);
    sqrtOneMinusSinSingleTest(-TWO_PI);
    for (int i = 0; i < 10000; i++) {
      final double a = java.lang.Math.random() * TWO_PI;
      sqrtOneMinusSinSingleTest(a);
      sqrtOneMinusSinSingleTest(-a);
    }
  }
  
  private static final void sqrtOneMinusSinSingleTest(final double a) {
    // System.out.println("aRad: " + a);
    double res = sqrtOfOneMinusSin(a);
    assertEquals(sqrt(1 - sin(a)), res, 1e-10);
    assertTrue(res >= 0);
    assertTrue(res <= SQRT2);
  }
  
  @Test
  public void sqrtOfOneMinusSinPCTest() {
    sqrtOneMinusSinPCSingleTest(0);
    sqrtOneMinusSinPCSingleTest(HALF_PI);
    sqrtOneMinusSinPCSingleTest(-HALF_PI);
    sqrtOneMinusSinPCSingleTest(-PI);
    for (int i = 0; i < 10000; i++) {
      final double a = -3 * HALF_PI + 4 * java.lang.Math.random() * HALF_PI;
      sqrtOneMinusSinPCSingleTest(a);
    }
  }
  
  private static final void sqrtOneMinusSinPCSingleTest(final double a) {
    // System.out.println("aRad: " + a);
    double res = sqrtOfOneMinusSinPC(a);
    assertEquals(sqrt(1 - sin(a)), res, 1e-10);
    assertTrue(res >= 0);
    assertTrue(res <= SQRT2);
  }
  
  
  public static void main(final String[] args) {
    final MathTest t = new MathTest();
    // t.oneMinusSinTest();
    // t.sqrtOfOneMinusSinTest();
    // t.sqrtOfOneMinusSinPCTest();
  }
}
