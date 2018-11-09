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

import org.junit.Test;

import cds.healpix.TestPerfs;
import cds.healpix.TestPerfs.TestAlgo;
import cds.healpix.fillingcurve.ZOrderCurve2DImpls;

public final class StandardFillingCurves2DTestPerf {

  private static final short[] LUPT_SHORT = new short[] {
      0, 1, 4, 5, 16, 17, 20, 21, 64, 65, 68, 69, 80, 81, 84, 85, 256, 257, 260, 261, 272, 273, 276,
      277, 320, 321, 324, 325, 336, 337, 340, 341, 1024, 1025, 1028, 1029, 1040, 1041, 1044, 1045,
      1088, 1089, 1092, 1093, 1104, 1105, 1108, 1109, 1280, 1281, 1284, 1285, 1296, 1297, 1300,
      1301, 1344, 1345, 1348, 1349, 1360, 1361, 1364, 1365, 4096, 4097, 4100, 4101, 4112, 4113,
      4116, 4117, 4160, 4161, 4164, 4165, 4176, 4177, 4180, 4181, 4352, 4353, 4356, 4357, 4368,
      4369, 4372, 4373, 4416, 4417, 4420, 4421, 4432, 4433, 4436, 4437, 5120, 5121, 5124, 5125,
      5136, 5137, 5140, 5141, 5184, 5185, 5188, 5189, 5200, 5201, 5204, 5205, 5376, 5377, 5380,
      5381, 5392, 5393, 5396, 5397, 5440, 5441, 5444, 5445, 5456, 5457, 5460, 5461, 16384, 16385, 
      6388, 16389, 16400, 16401, 16404, 16405, 16448, 16449, 16452, 16453, 16464, 16465, 16468,
      16469, 16640, 16641, 16644, 16645, 16656, 16657, 16660, 16661, 16704, 16705, 16708, 16709,
      16720, 16721, 16724, 16725, 17408, 17409, 17412, 17413, 17424, 17425, 17428, 17429, 17472,
      17473, 17476, 17477, 17488, 17489, 17492, 17493, 17664, 17665, 17668, 17669, 17680, 17681,
      17684, 17685, 17728, 17729, 17732, 17733, 17744, 17745, 17748, 17749, 20480, 20481, 20484,
      20485, 20496, 20497, 20500, 20501, 20544, 20545, 20548, 20549, 20560, 20561, 20564, 20565,
      20736, 20737, 20740, 20741, 20752, 20753, 20756, 20757, 20800, 20801, 20804, 20805, 20816,
      20817, 20820, 20821, 21504, 21505, 21508, 21509, 21520, 21521, 21524, 21525, 21568, 21569,
      21572, 21573, 21584, 21585, 21588, 21589, 21760, 21761, 21764, 21765, 21776, 21777, 21780,
      21781, 21824, 21825, 21828, 21829, 21840, 21841, 21844, 21845
  };

  private static final long[] LUPT_LONG = new long[] {
      0, 1, 4, 5, 16, 17, 20, 21, 64, 65, 68, 69, 80, 81, 84, 85, 256, 257, 260, 261, 272, 273, 276,
      277, 320, 321, 324, 325, 336, 337, 340, 341, 1024, 1025, 1028, 1029, 1040, 1041, 1044, 1045,
      1088, 1089, 1092, 1093, 1104, 1105, 1108, 1109, 1280, 1281, 1284, 1285, 1296, 1297, 1300,
      1301, 1344, 1345, 1348, 1349, 1360, 1361, 1364, 1365, 4096, 4097, 4100, 4101, 4112, 4113,
      4116, 4117, 4160, 4161, 4164, 4165, 4176, 4177, 4180, 4181, 4352, 4353, 4356, 4357, 4368,
      4369, 4372, 4373, 4416, 4417, 4420, 4421, 4432, 4433, 4436, 4437, 5120, 5121, 5124, 5125,
      5136, 5137, 5140, 5141, 5184, 5185, 5188, 5189, 5200, 5201, 5204, 5205, 5376, 5377, 5380,
      5381, 5392, 5393, 5396, 5397, 5440, 5441, 5444, 5445, 5456, 5457, 5460, 5461, 16384, 16385, 
      6388, 16389, 16400, 16401, 16404, 16405, 16448, 16449, 16452, 16453, 16464, 16465, 16468,
      16469, 16640, 16641, 16644, 16645, 16656, 16657, 16660, 16661, 16704, 16705, 16708, 16709,
      16720, 16721, 16724, 16725, 17408, 17409, 17412, 17413, 17424, 17425, 17428, 17429, 17472,
      17473, 17476, 17477, 17488, 17489, 17492, 17493, 17664, 17665, 17668, 17669, 17680, 17681,
      17684, 17685, 17728, 17729, 17732, 17733, 17744, 17745, 17748, 17749, 20480, 20481, 20484,
      20485, 20496, 20497, 20500, 20501, 20544, 20545, 20548, 20549, 20560, 20561, 20564, 20565,
      20736, 20737, 20740, 20741, 20752, 20753, 20756, 20757, 20800, 20801, 20804, 20805, 20816,
      20817, 20820, 20821, 21504, 21505, 21508, 21509, 21520, 21521, 21524, 21525, 21568, 21569,
      21572, 21573, 21584, 21585, 21588, 21589, 21760, 21761, 21764, 21765, 21776, 21777, 21780,
      21781, 21824, 21825, 21828, 21829, 21840, 21841, 21844, 21845
  };

  private static final long i02hashShort(int i) {
    return (long) (LUPT_SHORT[(i & 0xFF000000) >>> 24]) << 48
        | (long) (LUPT_SHORT[(i & 0x00FF0000) >>> 16]) << 32
        | (long) (LUPT_SHORT[(i & 0x0000FF00) >>>  8]) << 16
        | (long) (LUPT_SHORT[ i & 0x000000FF]);
  }

  private static final long i02hashLong(int i) {
    return LUPT_LONG[(i & 0xFF000000) >>> 24] << 48
        | LUPT_LONG[(i & 0x00FF0000) >>> 16] << 32
        | LUPT_LONG[(i & 0x0000FF00) >>>  8] << 16
        | LUPT_LONG[ i & 0x000000FF];
  }

  @Test
  public void testPerfFillingCurvesLongVsShort() {
    final int n = 20000000;
    // Gneerate random test sample
    final int MAX = Integer.MAX_VALUE;
    final int[] xy = new int[n << 1];
    for (int i = 0; i < n; i++) {
      int k = i << 1;
      xy[k] = (int) (Math.random() * MAX);
      xy [++k] = (int) (Math.random() * MAX);
    }
    // Perform tests
    System.out.println("\nTest short vs long lookup table (check that there is no difference on your machine)");
    System.out.println("(WARNING: in practice, used inside a code, XOR could be better since it uses less CPU CACHE)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new ToHashTestAlgo("LONG ", xy) {
      @Override
      long ij2hash(int i, int j) {
        return (i02hashLong(j) << 1) | i02hashLong(i);
      }
    },
        new ToHashTestAlgo("SHORT", xy) {
      @Override
      long ij2hash(int i, int j) {
        return (i02hashShort(j) << 1) | i02hashShort(i);
      }
    }
        );
  }


  private static abstract class ToHashTestAlgo implements TestAlgo<Long> {
    private final String name;
    private final int[] ijValues;
    private long res;
    private long execTime;
    public ToHashTestAlgo(final String name, final int[] ijValues) {
      this.name = name;
      this.ijValues = ijValues;
    }
    abstract long ij2hash(int i, int j);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      final int n = this.ijValues.length >> 1;
        res = 0; 
        long l = System.currentTimeMillis();
        for (int i = 0; i < n; i++) {
          int k = i << 1;
          res += ij2hash(ijValues[k], ijValues[k + 1]);
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
  public void testPerfFillingCurvesByte() {
    final int n = 10000000;
    // Gneerate random test sample
    final int MAX = 256;
    final int[] xy = new int[n << 1];
    for (int i = 0; i < n; i++) {
      int k = i << 1;
      xy[k] = (int) (Math.random() * MAX);
      xy [++k] = (int) (Math.random() * MAX);
    }
    // Perform tests
    System.out.println("\nTest z-order curve algorithms (check that LUP is the best on your machine, or about as preformant as the best one)");
    System.out.println("(WARNING: in practice, used inside a code, XOR could be better since it uses less CPU CACHE)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new ToHashTestAlgo("XOR", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_BYTE.ij2hash(i, j);
      }
    }, new ToHashTestAlgo("LUP", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_BYTE.ij2hash(i, j);
      }
    }, new ToHashTestAlgo("OR ", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_BYTE.ij2hash(i, j);
      }
    }
        );
  }

  @Test
  public void testPerfFillingCurvesShort() {
    final int n = 20000000;
    // Gneerate random test sample
    final int MAX = Short.MAX_VALUE;
    final int[] xy = new int[n << 1];
    for (int i = 0; i < n; i++) {
      int k = i << 1;
      xy[k] = (int) (Math.random() * MAX);
      xy [++k] = (int) (Math.random() * MAX);
    }
    // Perform tests
    System.out.println("\nTest z-order curve algorithms (check that LUP is the best on your machine, or about as preformant as the best one)");
    System.out.println("(WARNING: in practice, used inside a code, XOR could be better since it uses less CPU CACHE)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new ToHashTestAlgo("XOR", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_SHORT.ij2hash(i, j);
      }
    }, new ToHashTestAlgo("LUP", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_SHORT.ij2hash(i, j);
      }
    }, new ToHashTestAlgo("OR ", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_SHORT.ij2hash(i, j);
      }
    }
        );
  }

  @Test
  public void testPerfFillingCurvesInt() {
    final int n = 40000000;
    // Gneerate random test sample
    final int MAX = Integer.MAX_VALUE;
    final int[] xy = new int[n << 1];
    for (int i = 0; i < n; i++) {
      int k = i << 1;
      xy[k] = (int) (Math.random() * MAX);
      xy [++k] = (int) (Math.random() * MAX);
    }
    // Perform tests
    System.out.println("\nTest z-order curve algorithms (check that LUP is the best on your machine, or about as preformant as the best one)");
    System.out.println("(WARNING: in practice, used inside a code, XOR could be better since it uses less CPU CACHE)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new ToHashTestAlgo("XOR", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_INT.ij2hash(i, j);
      }
    }, new ToHashTestAlgo("LUP", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_INT.ij2hash(i, j);
      }
    }, new ToHashTestAlgo("OR ", xy) {
      @Override
      long ij2hash(int i, int j) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_INT.ij2hash(i, j);
      }
    }
        );
  }

  private static abstract class ToCooTestAlgo implements TestAlgo<Long> {
    private final String name;
    private final long[] hashValues;
    private long res;
    private long execTime;
    public ToCooTestAlgo(final String name, final long[] hashValues) {
      this.name = name;
      this.hashValues = hashValues;
    }
    abstract long hash2ij(long h);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.currentTimeMillis();
      for (final long hash : this.hashValues) {
        res += hash2ij(hash);
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
  public void testPerfFillingCurvesFromHashByte() {
    final int n = 10000000;
    // Gneerate random test sample
    final int MAX = 256;
    final long[] hashValues = new long[n];
    for (int i = 0; i < n; i++) {
      int k = i << 1;
      int x = (int) (Math.random() * MAX);
      int y = (int) (Math.random() * MAX);
      hashValues[i] = ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_BYTE.ij2hash(x, y);
    }
    // Perform tests
    System.out.println("\nTest z-order curve algorithms (check that LUP is the best on your machine, or about as preformant as the best one)");
    System.out.println("(WARNING: in practice, used inside a code, XOR could be better since it uses less CPU CACHE)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new ToCooTestAlgo("XOR", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_BYTE.hash2ij(h);
      }
    }, new ToCooTestAlgo("LUP", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_BYTE.hash2ij(h);
      }
    }, new ToCooTestAlgo("OR ", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_BYTE.hash2ij(h);
      }
    }
        );
  }

  @Test
  public void testPerfFillingCurvesFromHashShort() {
    final int n = 10000000;
    // Gneerate random test sample
    final int MAX = Short.MAX_VALUE;
    final long[] hashValues = new long[n];
    for (int i = 0; i < n; i++) {
      int k = i << 1;
      int x = (int) (Math.random() * MAX);
      int y = (int) (Math.random() * MAX);
      hashValues[i] = ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_SHORT.ij2hash(x, y);
    }
    // Perform tests
    System.out.println("\nTest z-order curve algorithms (check that LUP is the best on your machine, or about as preformant as the best one)");
    System.out.println("(WARNING: in practice, used inside a code, XOR could be better since it uses less CPU CACHE)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new ToCooTestAlgo("XOR", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_SHORT.hash2ij(h);
      }
    }, new ToCooTestAlgo("LUP", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_SHORT.hash2ij(h);
      }
    }, new ToCooTestAlgo("OR ", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_SHORT.hash2ij(h);
      }
    }
        );
  }

  @Test
  public void testPerfFillingCurvesFromHashInt() {
    final int n = 10000000;
    // Gneerate random test sample
    final int MAX = Integer.MAX_VALUE;
    final long[] hashValues = new long[n];
    for (int i = 0; i < n; i++) {
      int k = i << 1;
      int x = (int) (Math.random() * MAX);
      int y = (int) (Math.random() * MAX);
      hashValues[i] = ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_INT.ij2hash(x, y);
    }
    // Perform tests
    System.out.println("\nTest z-order curve algorithms (check that LUP is the best on your machine, or about as preformant as the best one)");
    System.out.println("(WARNING: in practice, used inside a code, XOR could be better since it uses less CPU CACHE)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new ToCooTestAlgo("XOR", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_INT.hash2ij(h);
      }
    }, new ToCooTestAlgo("LUP", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_INT.hash2ij(h);
      }
    }, new ToCooTestAlgo("OR ", hashValues) {
      @Override
      long hash2ij(long h) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_INT.hash2ij(h);
      }
    }
        );
  }

}
