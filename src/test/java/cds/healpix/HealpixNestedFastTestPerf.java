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

import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.toBits;

import java.util.Random;

import org.junit.Test;

import cds.healpix.TestPerfs.TestAlgo;

public final class HealpixNestedFastTestPerf {

  private static abstract class TimeHalfNiseTestAlgo implements TestAlgo<Long> {
    private final String name;
    private final double[] values;
    private long res;
    private long execTime;
    public TimeHalfNiseTestAlgo(final String name, final double[] values) {
      this.name = name;
      this.values = values;
    }
    abstract long timeHalfNside(double v);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.currentTimeMillis();
      for (final double v : this.values) {
        res += timeHalfNside(v);
      }
      execTime = System.currentTimeMillis() - l;
    }
    @Override
    public Long getResult() {
      return res;
    }
    @Override
    public long getExecutionDuration() {
      return execTime;
    }
  }
  
  @Test
  public void timeHalfNsideTestPerf() {
    final int depth = 12;
    final int nTests = 20000000;
    // final long nTest = Healpix.nHash(depth);
    final long longHalfNside = (depth - 1L) << 52;
    final double halfNside = (1 << depth) >> 1;
    final double[] sample = new double[nTests];
    for (int i = 0 ; i < nTests; i++) {
      sample[i] = Math.random() * 4;
    }
    System.out.println("Test timeHalfNside algorithms (check that IEEEv2 is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(3, new TimeHalfNiseTestAlgo("IEEEv1", sample) {
      @Override
      long timeHalfNside(double v) {
        long bits = toBits(v);
        long shift = 1076 - (depth + ((bits & 0x7FF0000000000000L) >>> 52));
        return shift > 52 || shift < 0 ? 0 : ((bits & 0x000FFFFFFFFFFFFFL) | 0x0010000000000000L) >>> shift;
      }
    }, new TimeHalfNiseTestAlgo("IEEEv2", sample) {
      @Override
      long timeHalfNside(double v) {
        return (long) fromBits(toBits(v) + longHalfNside);
      }
    }, new TimeHalfNiseTestAlgo("CLASSIC", sample) {
      @Override
      long timeHalfNside(double v) {
        return (long) (v *  halfNside);
      }
    }
    );
    System.out.println("Test timeHalfNside algorithms (check that IEEEv2 is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(3, new TimeHalfNiseTestAlgo("CLASSIC", sample) {
      @Override
      long timeHalfNside(double v) {
        return (long) (v *  halfNside);
      }
    }, new TimeHalfNiseTestAlgo("IEEEv1", sample) {
      @Override
      long timeHalfNside(double v) {
        long bits = toBits(v);
        long shift = 1076 - (depth + ((bits & 0x7FF0000000000000L) >>> 52));
        return shift > 52 || shift < 0 ? 0 : ((bits & 0x000FFFFFFFFFFFFFL) | 0x0010000000000000L) >>> shift;
      }
    }, new TimeHalfNiseTestAlgo("IEEEv2", sample) {
      @Override
      long timeHalfNside(double v) {
        return (long) fromBits(toBits(v) + longHalfNside);
      }
    }
    );
    System.out.println("Test timeHalfNside algorithms (check that IEEEv2 is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(3, new TimeHalfNiseTestAlgo("IEEEv2", sample) {
      @Override
      long timeHalfNside(double v) {
        return (long) fromBits(toBits(v) + longHalfNside);
      }
    }, new TimeHalfNiseTestAlgo("CLASSIC", sample) {
      @Override
      long timeHalfNside(double v) {
        return (long) (v *  halfNside);
      }
    }, new TimeHalfNiseTestAlgo("IEEEv1", sample) {
      @Override
      long timeHalfNside(double v) {
        long bits = toBits(v);
        long shift = 1076 - (depth + ((bits & 0x7FF0000000000000L) >>> 52));
        return shift > 52 || shift < 0 ? 0 : ((bits & 0x000FFFFFFFFFFFFFL) | 0x0010000000000000L) >>> shift;
      }
    }
    );
  }

  
}
