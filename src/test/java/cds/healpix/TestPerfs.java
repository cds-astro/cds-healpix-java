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

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public final class TestPerfs {

  @Test
  public void dummy() {
    assertTrue(true); // Let ant be happy testing this class
  }
  
  public static interface TestAlgo<E extends Comparable<? super E>> {
    String algoName();
    void exec();
    long getExecutionDuration();
    E getResult();
  }

  // @SafeVarargs
  public static <T extends Comparable<? super T>> void test(int n, TestAlgo<T>... algos) {
    for (int i = 0; i < n; i++) {
      performTests(algos);
      compareResults(algos);
      printResults(algos);
    }
  }

  private static void performTests(final TestAlgo<?>[] algos) {
    for (final TestAlgo<?> algo : algos) {
      algo.exec();
    }
  }
  
  private static <T extends Comparable<? super T>> void compareResults(final TestAlgo<T>[] algos) {
    T ref = algos[0].getResult();
    for (int i = 1; i < algos.length; i++) {
      if (ref.compareTo(algos[i].getResult()) != 0) {
        throw new Error(ref + " != " + algos[i].getResult());
      }
    }
  }

  private static void printResults(final TestAlgo<?>[] algos) {
    final TestAlgo<?> best = findQuickest(algos);
    System.out.println("---------");
    System.out.println("Best algo: " + best.algoName() + "; exec in: " + best.getExecutionDuration() + " ms");
    for (final TestAlgo<?> algo : algos) {
      int percent = (int)(100 * algo.getExecutionDuration() / (double) best.getExecutionDuration());
      System.out.println(" - " + algo.algoName() + " exec in " + algo.getExecutionDuration()
          + " ms => " + percent + " % of best algo");
    }
  }

  private static TestAlgo<?> findQuickest(final TestAlgo<?>[] algos) {
    int iQuickest = 0;
    long minExecTime = algos[0].getExecutionDuration();
    for (int i = 1; i < algos.length; i++) {
      long t = algos[i].getExecutionDuration();
      if (t < minExecTime) {
        iQuickest = i;
        minExecTime = t;
      }
    }
    return algos[iQuickest];
  }
}
