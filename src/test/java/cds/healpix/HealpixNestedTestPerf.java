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

import org.junit.Test;

import cds.healpix.TestPerfs.TestAlgo;

public class HealpixNestedTestPerf {

  private static abstract class BasePixelTestAlgo implements TestAlgo<Integer> {
    private final String name;
    private final int[][] ijValues;
    private int res;
    private long execTime;
    public BasePixelTestAlgo(final String name, final int[][] ijValues) {
      this.name = name;
      this.ijValues = ijValues;
    }
    abstract int basePixel(int i, int j);
    @Override
    public String algoName() {
      return this.name;
    }
    @Override
    public void exec() {
      res = 0; 
      long l = System.nanoTime();
      for (final int[] ij : this.ijValues) {
        res += basePixel(ij[0], ij[1]);
      }
      execTime = System.nanoTime() - l;
    }
    @Override
    public Integer getResult() {
      return res;
    }
    @Override
    public long getExecutionDuration() {
      return (long) (execTime / 1000000.0);
    }
  }
  @Test
  public void basePixelTestPerf() {
    final int[][] reverse = new int[][]{{1, 0}, {2, -1}, {-1, 2}, {0, 1},
                                        {0, 0}, {1, -1}, {2, -2}, {-1, 1},
                                        {0, -1}, {1, -2}, {-2, 1}, {-1, 0}};
    final int[] LOOKUP1 = new int[]{4, 3, -1, 8, 0, -1, 9, 5, -1, 10, 6, 1, 11, 7, 2, -1};
    final int[][] LOOKUP2 = new int[][]{{4, 3, -1, 8}, {0, -1, 9, 5}, {-1, 10, 6, 1}, {11, 7, 2, -1}};
    final int n = 10000000;
    int[][] coos = new int[n][2];
    // Prepare random nside values
    for (int i = 0; i < n; i++) {
      int basePix = (int) (Math.random() * 12);
      coos[i] = reverse[basePix];
    }
    // Run tests
    System.out.println("\nTest basePixel algorithms (check that NOBRANCH is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new BasePixelTestAlgo("NOBRANCH", coos)  {
      @Override
      int basePixel(int i, int j) {
        int offset = (i - (i << 1) - j);
        i -= (offset >>> 63);
        offset++;
        offset <<= 2;
        i &= 3;
        return i + offset;
      }
    }, new BasePixelTestAlgo("LOOKUPT1", coos)  {
      @Override
      int basePixel(int i, int j) {
        return LOOKUP1[((i & 3) << 2) + (j & 3)];
      }
    }, new BasePixelTestAlgo("LOOKUPT2", coos)  {
      @Override
      int basePixel(int i, int j) {
        return LOOKUP2[i & 3][j & 3];
      }
    }
    );
  }
  
  private static final byte[][] D0H_LOOKUP_MATRIX = new byte[][]{        
    {-1, -1, -1,  8,  4}, //   ----> y-axis 
    {-1, -1,  9,  5,  0}, //  |
    {-1, 10,  6,  1, -1}, //  |
    {11,  7,  2, -1, -1}, //  v
    { 4,  3, -1, -1, -1}  // x-axis
  };
  
  private static final byte[] D0H_LOOKUP_ARRAY = new byte[]{        
    -1, -1, -1,  8,  4, //   ----> y-axis 
    -1, -1,  9,  5,  0, //  |
    -1, 10,  6,  1, -1, //  |
    11,  7,  2, -1, -1, //  v
     4,  3, -1, -1, -1  // x-axis
  };

  @Test
  public void basePixelTestPerf2() {
    final int n = 50000000;
    int[][] coos = new int[n][2];
    // Prepare random nside values
    int ai = 0;
    while (ai < n) {
      int i = (int) (Math.random() * 5);
      int j = (int) (Math.random() * 5);
      if (i < 5 && j < 5 && D0H_LOOKUP_MATRIX[i][j] != -1) {
        coos[ai]= new int[]{i, j};
        ai++;
      }
    }
    // Run tests
    System.out.println("\nTest basePixel algorithms (check that LOOKUP_MATRIX is the best on your machine, or about as preformant as the best one)");
    System.out.println("(Rmk: there is a same overhead for each algo due to abstraction in testing)");
    TestPerfs.test(4, new BasePixelTestAlgo("NOBRANCH", coos)  {
      @Override
      int basePixel(int i, int j) {
        j = 4 - (j + i);
        return ((i - (j >>> 63)) & 3) + (++j << 2);
      }
    }, new BasePixelTestAlgo("LOOKUP_MATRIX", coos)  {
      @Override
      int basePixel(int i, int j) {
        return D0H_LOOKUP_MATRIX[i][j];
      }
    }, new BasePixelTestAlgo("LOOKUP_ARRAY", coos)  {
      @Override
      int basePixel(int i, int j) {
        return D0H_LOOKUP_ARRAY[i * 5 + j];
      }
    }, new BasePixelTestAlgo("LOOKUP_ARRAY v2", coos)  {
      @Override
      int basePixel(int i, int j) {
        return D0H_LOOKUP_ARRAY[i + j + (i << 2)];
      }
    }
    );
  }
}
