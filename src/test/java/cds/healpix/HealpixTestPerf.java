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

import static cds.healpix.common.math.HackersDelight.BUT_SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.nlz;
import static cds.healpix.common.math.HackersDelight.toBits;
import static cds.healpix.common.math.Math.ONE_OVER_PI;
import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.TWO_OVER_PI;
import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.log;

import java.util.Arrays;

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.TestPerfs.TestAlgo;


public final class HealpixTestPerf {

  private static abstract class AbsTestAlgo implements TestAlgo<Double> {
    private final String name;
    private final double[] values;
    private double res;
    private long execTime;
    public AbsTestAlgo(final String name, final double[] values) {
      this.name = name;
      this.values = values;
    }
    abstract double abs(double val);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.nanoTime();
      for (double value : this.values) {
        res += abs(value);
      }
      execTime = System.nanoTime() - l;
    }
    @Override
    public Double getResult() {
      return res;
    }
    @Override
    public long getExecutionDuration() {
      return execTime / 1000000;
    }
  }

  @Test
  public void absTestPerf() {
    final int n = 25000000;
    double[] longitudes = new double[n];
    // Prepare random nside values
    for (int i = 0; i < n; i++) {
      longitudes[i] = (Math.random() - 0.5);
      /*if (i < 100) {
        System.out.println(longitudes[i]);
      }*/
    }
    // Run tests
    System.out.println("Test abs algorithms (check that BITS1 is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new AbsTestAlgo("BITS1", longitudes) {
      @Override
      double abs(final double val) {
        return fromBits(toBits(val) & BUT_SIGN_BIT_MASK_L);
      }
    },
        new AbsTestAlgo("BITS2", longitudes) {
      @Override
      double abs(final double val) {
        return fromBits((toBits(val) << 1) >>> 1);
      }
    },
        new AbsTestAlgo("BITS3", longitudes) {
      @Override
      double abs(final double val) {
        long l = (toBits(val) & SIGN_BIT_MASK_L) >> 62;
        return ++l * val;
      }
    },
        new AbsTestAlgo("TEST ", longitudes) {
      @Override
      double abs(final double val) {
        return val <= 0d ? 0d - val : val;
      }
    },
        new AbsTestAlgo("MATH ", longitudes) {
      @Override
      double abs(final double val) {
        return java.lang.Math.abs(val);
      }
    }/*,
        new AbsTestAlgo("FMATH", longitudes) {
      @Override
      double abs(final double val) {
        return org.apache.commons.math3.util.FastMath.abs(val);
      }
    }*/
        );
  }


  private static abstract class DepthTestAlgo implements TestAlgo<Long> {
    private final String name;
    private final int[] nsides;
    private long res;
    private long execTime;
    public DepthTestAlgo(final String name, final int nsides[]) {
      this.name = name;
      this.nsides = nsides;
    }
    abstract int depth(int nside);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.nanoTime();
      for (int nside : this.nsides) {
        res += depth(nside);
      }
      execTime = System.nanoTime() - l;
    }
    @Override
    public Long getResult() {
      return res;
    }
    @Override
    public long getExecutionDuration() {
      return execTime / 1000000;
    }
  }

  @Test
  public void depthTestPerf() {
    final int n = 20000000;
    int[] nsides = new int[n];
    // Prepare random nside values
    for (int i = 0; i < n; i++) {
      nsides[i] = Healpix.nside((int) (Math.random() * 29));
    }
    // Run tests
    System.out.println("Test depth algorithms (check that TRAIL is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4,
        new DepthTestAlgo("WHILE", nsides) {
      @Override
      int depth(int nside) {
        int i = 0;
        while ((nside >> (++i)) > 0);
        return --i;
      }
    }, new DepthTestAlgo("MATH ", nsides) {
      @Override
      int depth(int nside) {
        return (int) (Math.log(nside) / Math.log(2));
      }
    }, new DepthTestAlgo("TRAIL", nsides) {
      @Override
      int depth(int nside) {
        return Integer.numberOfTrailingZeros(nside);
      }
    }, new DepthTestAlgo("FASTM", nsides) {
      @Override
      int depth(int nside) {
        return (int) (log(nside) / log(2));
      }
    }, new DepthTestAlgo("NZL  ", nsides) {
      @Override
      int depth(int nside) {
        return 31 - nlz(nside);
      }
    }
        );
  } 


  private static abstract class TestTestAlgo implements TestAlgo<Long> {
    private final String name;
    private final double[] values;
    private long nPos;
    private long execTime;
    public TestTestAlgo(final String name, final double[] values) {
      this.name = name;
      this.values = values;
    }
    abstract boolean isInEquatorialRegion(double val);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      nPos = 0; 
      long l = System.nanoTime();
      for (double value : this.values) {
        if (isInEquatorialRegion(value)) {
          nPos++;
        }
      }
      execTime = System.nanoTime() - l;
    }
    @Override
    public Long getResult() {
      return nPos;
    }
    @Override
    public long getExecutionDuration() {
      return execTime / 1000000;
    }
  }
  @Test
  public void testTestPerf() {
    final int n = 40000000;
    double[] longitudes = new double[n];
    // Prepare random nside values
    for (int i = 0; i < n; i++) {
      longitudes[i] = (Math.random() - 0.5);
    }
    // Run tests
    System.out.println("Test isInEquatorialRegion algorithms (check that ABS1 is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    final double limit = 0.25;
    TestPerfs.test(4, new TestTestAlgo("ABS1", longitudes) {
      @Override
      boolean isInEquatorialRegion(final double val) {
        return fromBits(toBits(val) & BUT_SIGN_BIT_MASK_L) <= limit;
      }
    }, new TestTestAlgo("TEST1", longitudes) {
      @Override
      boolean isInEquatorialRegion(final double val) {
        return -limit <= val && val <= limit;
      }
    },new TestTestAlgo("ABS2 ", longitudes) {
      @Override
      boolean isInEquatorialRegion(final double val) {
        return fromBits((toBits(val) << 1) >>> 1) <= limit;
      }
    }, new TestTestAlgo("TEST2", longitudes) {
      @Override
      boolean isInEquatorialRegion(final double val) {
        return !(val < -limit || val > limit);
      }
    }
        );
  }

  private static abstract class NormalizeFromMinusPiToPiTestAlgo implements TestAlgo<Double> {
    private final String name;
    private final double[] values;
    private double res;
    private long execTime;
    public NormalizeFromMinusPiToPiTestAlgo(final String name, final double[] values) {
      this.name = name;
      this.values = values;
    }
    abstract double normalizeLonFromMinusPiToPi(double val);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.nanoTime();
      for (double value : this.values) {
        res += normalizeLonFromMinusPiToPi(value);
      }
      execTime = System.nanoTime() - l;
    }
    @Override
    public Double getResult() {
      return res;
    }
    @Override
    public long getExecutionDuration() {
      return execTime / 1000000;
    }
  }
  @Test
  public void normalizeLonFromMinusPiToPiTestPerf() {
    final int n = 80000000;
    double[] longitudes = new double[n];
    // Prepare random nside values
    for (int i = 0; i < n; i++) {
      longitudes[i] = Math.random() * TWO_PI;
    }
    Arrays.sort(longitudes);
    // Run tests
    System.out.println("Test normalizeLonFromMinusPiToPi algorithms (check that COMP is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    final double[] OFF = new double[]{0, TWO_PI};
    final long TWO_PI_BITS = toBits(TWO_PI);
    TestPerfs.test(4, new NormalizeFromMinusPiToPiTestAlgo("TEST1", longitudes) {
      @Override
      double normalizeLonFromMinusPiToPi(double lon) {
        return lon > PI ? lon - TWO_PI : lon;
      }
    }, new NormalizeFromMinusPiToPiTestAlgo("COMP ", longitudes) {
      @Override
      double normalizeLonFromMinusPiToPi(double lon) {
        int half = (int) (lon * ONE_OVER_PI);
        lon -= half * TWO_PI;
        return lon;
      }
    }, new NormalizeFromMinusPiToPiTestAlgo("BITS ", longitudes) {
      @Override
      double normalizeLonFromMinusPiToPi(double lon) {
        long half = (long) (lon * ONE_OVER_PI); // 0 ou 1
        lon -= fromBits(TWO_PI_BITS & ~(--half));
        return lon;
      }
    }, new NormalizeFromMinusPiToPiTestAlgo("MODULO ", longitudes) {
      @Override
      double normalizeLonFromMinusPiToPi(double lon) {
        int n = (int) (lon % PI);
        long half = (long) (lon * ONE_OVER_PI); // 0 ou 1
        lon -= fromBits(TWO_PI_BITS & ~(--half));
        return lon;
      }
    }, new NormalizeFromMinusPiToPiTestAlgo("TEST2", longitudes) {
      @Override
      double normalizeLonFromMinusPiToPi(double lon) {
        return lon > PI ? lon - TWO_PI : lon;
      }
    }, new NormalizeFromMinusPiToPiTestAlgo("LUPT ", longitudes) {
      @Override
      double normalizeLonFromMinusPiToPi(double lon) {
        return lon -  OFF[(int) (lon * ONE_OVER_PI)];
      }
    }
        );
  }


  private static abstract class NormalizeFrom0ToTwicePiTestAlgo implements TestAlgo<Double> {
    private final String name;
    private final double[] values;
    private double res;
    private long execTime;
    public NormalizeFrom0ToTwicePiTestAlgo(final String name, final double[] values) {
      this.name = name;
      this.values = values;
    }
    abstract double normalizeLonFrom0ToTwicePi(double val);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.nanoTime();
      for (double value : this.values) {
        res += normalizeLonFrom0ToTwicePi(value);
      }
      execTime = System.nanoTime() - l;
    }
    @Override
    public Double getResult() {
      return res;
    }
    @Override
    public long getExecutionDuration() {
      return execTime / 1000000;
    }
  }
  @Test
  public void normalizeLonFrom0ToTwicePiTestPerf() {
    final int n = 40000000;
    double[] longitudes = new double[n];
    // Prepare random nside values
    for (int i = 0; i < n; i++) {
      longitudes[i] = Math.random() * TWO_PI - PI;
    }
    // Run tests
    System.out.println("Test normalizeLonFrom0ToTwicePi algorithms (check that BITS2 is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    final double[] OFF = new double[]{0, TWO_PI};
    final long TWO_PI_BITS = toBits(TWO_PI);
    TestPerfs.test(4, new NormalizeFrom0ToTwicePiTestAlgo("BITS1", longitudes) {
      @Override
      double normalizeLonFrom0ToTwicePi(double lon) {
        long sign = toBits(lon);
        sign >>>= 63;
        lon += sign * TWO_PI;
        return lon;
      }
    }, new NormalizeFrom0ToTwicePiTestAlgo("BITS2", longitudes) {
      @Override
      double normalizeLonFrom0ToTwicePi(double lon) {
        long sign = toBits(lon);
        sign >>= 63;
    lon += fromBits(TWO_PI_BITS & sign);
      return lon;
      }
    }, new NormalizeFrom0ToTwicePiTestAlgo("TEST ", longitudes) {
      @Override
      double normalizeLonFrom0ToTwicePi(double lon) {
        return lon < 0 ? lon + TWO_PI : lon;
      }
    }, new NormalizeFrom0ToTwicePiTestAlgo("LUPT ", longitudes) {
      @Override
      double normalizeLonFrom0ToTwicePi(double lon) {
        return lon +  OFF[(int) (toBits(lon) >>> 63)];
      }
    }
        );
  }


  private static abstract class ComputeOffsetConstantTestAlgo implements TestAlgo<Double> {
    private final String name;
    private final double[] values;
    private double res;
    private long execTime;
    public ComputeOffsetConstantTestAlgo(final String name, final double[] values) {
      this.name = name;
      this.values = values;
    }
    abstract double computeOffsetConstant(double val);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.nanoTime();
      for (double value : this.values) {
        res += computeOffsetConstant(value);
      }
      execTime = System.nanoTime() - l;
    }
    @Override
    public Double getResult() {
      return res;
    }
    @Override
    public long getExecutionDuration() {
      return execTime / 1000000;
    }
  }
  @Test
  public void computeOffsetConstantTestPerf() {
    final int n = 80000000;
    double[] longitudes = new double[n];
    // Prepare random nside values
    for (int i = 0; i < n; i++) {
      longitudes[i] = Math.random() * TWO_PI;
    }
    // Run tests
    System.out.println("Test normalizeLonFrom0ToTwicePi algorithms (check that BITS2 is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    final double[] OFFSET_LUPT = new double[]{
        0.78539816339744830961, // 0 * pi / 2 + pi / 4 =   pi / 4 = 0.25 pi
        2.35619449019234492883, // 1 * pi / 2 + pi / 4 = 3 pi / 4 = 0.75 pi
        3.92699081698724154805, // 2 * pi / 2 + pi / 4 = 5 pi / 4 = 1.25 pi
        5.49778714378213816727  // 3 * pi / 2 + pi / 4 = 7 pi / 4 = 1.75 pi
    };
    TestPerfs.test(4, new ComputeOffsetConstantTestAlgo("LUPT", longitudes) {
      @Override
      double computeOffsetConstant(double lon) {              
        final int k = (int) (lon * TWO_OVER_PI);
        return OFFSET_LUPT[k];
      }
    }, new ComputeOffsetConstantTestAlgo("V1", longitudes) {
      @Override
      double computeOffsetConstant(double lon) {
        final int k = (int) (lon * TWO_OVER_PI);
        return 0.25 * PI * (2 * k + 1);
      }
    }, new ComputeOffsetConstantTestAlgo("V2 ", longitudes) {
      @Override
      double computeOffsetConstant(double lon) {
        final int k = (int) (lon * TWO_OVER_PI);
        return (PI / 4) * (2 * k + 1);
      }
    }, new ComputeOffsetConstantTestAlgo("V3 ", longitudes) {
      @Override
      double computeOffsetConstant(double lon) {
        int k = (int) (lon * TWO_OVER_PI);
        k <<= 1;
        ++k;
        return PI_OVER_FOUR * k;
      }
    }
        );
  }

}
