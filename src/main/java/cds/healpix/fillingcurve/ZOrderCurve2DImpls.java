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

package cds.healpix.fillingcurve;

/**
 * Utility class containing several implementations of the z-order curve in two dimension.
 * All implementations implement the {@link FillingCurve2D} interface.
 *
 * @author F.-X. Pineau
 *
 */
public final class ZOrderCurve2DImpls {

  /**
   * Private constructor in a final class prevents from instantiation outside this class.
   */
  private ZOrderCurve2DImpls() { }

  /**
   * Implementation doing nothing (all methods return 0).
   */
  public static final FillingCurve2D EMPTY = new FillingCurve2D() {
    @Override
    public long xy2hash(double x, double y) {
      return 0;
    }
    @Override
    public long ij2hash(int i, int j) {
      return 0;
    }
    @Override
    public long i02hash(int i) {
      return 0;
    }
    @Override
    public long hash2ij(long hash) {
      return 0;
    }
    @Override
    public long hash2i0(long hash) {
      return 0;
    }
    @Override
    public int ij2i(long ij) {
      return 0;
    }
    @Override
    public int ij2j(long ij) {
      return 0;
    }
  };
  
  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on the bitwise OR operator to interleave the bits of
   * the discretized 2d-coordinates. We assume that each discritized coordinates is coded on maximum
   * 8 bits (BYTE).
   * The algorithm is a slightly adapted version of the outer perfect shuffle define
   * p. 106 of "Hacker's Delight" (Henry S. Warren, Jr), slightly modified.
   */
  public static final FillingCurve2D ZOC_VMSB_OR_BYTE = new FillingCurve2D() {

    @Override
    public final long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public final long ij2hash(int i, int j) {
      // Verify for each coordinate that no bit is set after the 8th LSB
      assert (i == -1 && j == -1) || ((i & 0xFFFFFF00) == 0 && (j & 0xFFFFFF00) == 0);
      i |= (j << 8);
      i = (0x000000F0 & i) <<  4
          | (0x00000F00 & i) >>> 4
          | (0x0000F00F & i);
      i = (0x00000C0C & i) <<  2
          | (0x00003030 & i) >>> 2
          | (0x0000C3C3 & i);
      i = (0x00002222 & i) <<  1
          | (0x00004444 & i) >>> 1
          | (0x00009999 & i);
      return i;
    }

    @Override
    public final long i02hash(int i) {
      // Verify that the coordinate has no bit set after the 8th LSB
      assert (i & 0xFFFFFF00) == 0;
      // i = ((0x000000F0 & i) << 4) | (0x0000000F & i);
      i = ((i << 4) | i) & 0x00000F0F;
      i = ((i << 2) | i) & 0x00003333;
      i = ((i << 1) | i) & 0x00005555;
      return i;
    }

    @Override
    public final long hash2ij(long hash) {
      // Verify no bit set after the 16th LSB
      assert (0xFFFFFFFFFFFF0000L & hash) == 0;
      int h = (int) hash;
      h = (0x22222222 & h) <<  1
          | (0x44444444 & h) >>> 1
          | (0x99999999 & h);
      h = (0x0C0C0C0C & h) <<  2
          | (0x30303030 & h) >>> 2
          | (0xC3C3C3C3 & h);
      h = (0x00F000F0 & h) <<  4
          | (0x0F000F00 & h) >>> 4
          | (0xF00FF00F & h);
      return h;
    }

    @Override
    public final long hash2i0(long hash) {
      // Verify all bit of j are set to 0
      assert (0xFFFFFFFFFFFF3333L & hash) == 0;
      int h = (int) hash;
      h = ((h >> 1) | h) & 0x00003333;
      h = ((h >> 2) | h) & 0x00000F0F;
      h = ((h >> 4) | h) & 0x000000FF;
      return h;
    }

    @Override
    public final int ij2i(final long ij) {
      return (((int) ij) & 0x0000FF);
    }

    @Override
    public final int ij2j(final long ij) {
      return (((int) ij) >>> 8);
    }
  };

  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on the bitwise OR operator to interleave the bits of
   * the discretized 2d-coordinates. We assume that each discritized coordinates is coded on maximum
   * 16 bits (SHORT).
   * The algorithm is a slightly adapted version of the outer perfect shuffle define
   * p. 106 of "Hacker's Delight" (Henry S. Warren, Jr), slightly modified.
   */
  public static final FillingCurve2D ZOC_VMSB_OR_SHORT = new FillingCurve2D() {

    @Override
    public final long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public final long ij2hash(int i, int j) {
      // Verify for each coordinate that no bit is set after the 16th LSB
      assert (i == -1 && j == -1) || ((i & 0xFFFF0000) == 0 && (j & 0xFFFF0000) == 0);
      long h = ((long) i) | ((long) j << 16);
      // It is important to work on long because of sign problems when short coded in j is negative
      // i.e. when the 16^th it of of j = 1
      h = (0x000000000000FF00L & h) <<  8
          | (0x0000000000FF0000L & h) >>> 8
          | (0x00000000FF0000FFL & h);
      h = (0x0000000000F000F0L & h) <<  4
          | (0x000000000F000F00L & h) >>> 4
          | (0x00000000F00FF00FL & h);
      h = (0x000000000C0C0C0CL & h) <<  2
          | (0x0000000030303030L & h) >>> 2
          | (0x00000000C3C3C3C3L & h);
      h = (0x0000000022222222L & h) <<  1
          | (0x0000000044444444L & h) >>> 1
          | (0x0000000099999999L & h);
      return h;
    }

    @Override
    public final long i02hash(int i) {
      // Verify that the coordinate has no bit set after the 16th LSB
      assert (i & 0xFFFF0000) == 0;
      // i = ((0x0000FF00 & i) << 8) | (0x000000FF & i);
      i = ((i << 8) | i) & 0x00FF00FF;
      i = ((i << 4) | i) & 0x0F0F0F0F;
      i = ((i << 2) | i) & 0x33333333;
      i = ((i << 1) | i) & 0x55555555;
      return i;
    }

    @Override
    public final long hash2ij(long hash) {
      // Verify no bit set after the 32th LSB
      assert (0xFFFFFFFF00000000L & hash) == 0;
      int h = (int) hash;
      h = (0x22222222 & h) <<  1
          | (0x44444444 & h) >>> 1
          | (0x99999999 & h);
      h = (0x0C0C0C0C & h) <<  2
          | (0x30303030 & h) >>> 2
          | (0xC3C3C3C3 & h);
      h = (0x00F000F0 & h) <<  4
          | (0x0F000F00 & h) >>> 4
          | (0xF00FF00F & h);
      h = (0x0000FF00 & h) <<  8
          | (0x00FF0000 & h) >>> 8
          | (0xFF0000FF & h);
      return h;
    }

    @Override
    public final long hash2i0(long hash) {
      // Verify all bit of j are set to 0
      assert (0xFFFFFFFF33333333L & hash) == 0;
      int h = (int) hash;
      h = ((h >> 1) | h) & 0x33333333;
      h = ((h >> 2) | h) & 0x0F0F0F0F;
      h = ((h >> 4) | h) & 0x00FF00FF;
      h = ((h >> 8) | h) & 0x0000FFFF;
      return h;
    }

    @Override
    public final int ij2i(final long ij) {
      return (((int) ij) & 0x0000FFFF);
    }

    @Override
    public final int ij2j(final long ij) {
      return (((int) ij) >>> 16);
    }
  };

  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on the bitwise OR operator to interleave the bits of
   * the discretized 2d-coordinates. We assume that each discritized coordinates is coded on maximum
   * 32 bits.
   * The algorithm is a slightly adapted and extended version of the outer perfect shuffle define
   * p. 106 of "Hacker's Delight" (Henry S. Warren, Jr).
   */
  public static final FillingCurve2D ZOC_VMSB_OR_INT = new FillingCurve2D() {

    @Override
    public final long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public final long ij2hash(int i, int j) {
      long h = ((long) j ) << 32;
      h |= (long) i;
      h = (0x00000000FFFF0000L & h) <<  16
          | (0x0000FFFF00000000L & h) >>> 16
          | (0xFFFF00000000FFFFL & h);
      h = (0x0000FF000000FF00L & h) <<  8
          | (0x00FF000000FF0000L & h) >>> 8
          | (0xFF0000FFFF0000FFL & h);
      h = (0x00F000F000F000F0L & h) <<  4
          | (0x0F000F000F000F00L & h) >>> 4
          | (0xF00FF00FF00FF00FL & h);
      h = (0x0C0C0C0C0C0C0C0CL & h) <<  2
          | (0x3030303030303030L & h) >>> 2
          | (0xC3C3C3C3C3C3C3C3L & h);
      h = (0x2222222222222222L & h) <<  1
          | (0x4444444444444444L & h) >>> 1
          | (0x9999999999999999L & h);
      return h;
    }

    @Override
    public final long i02hash(int i) {
      long h = (long) i;
      /// h = ((0x00000000FFFF0000L & h) << 16) | (0x000000000000FFFFL & h);
      h = ((h << 16) | h) & 0x0000FFFF0000FFFFL;
      h = ((h <<  8) | h) & 0x00FF00FF00FF00FFL;
      h = ((h <<  4) | h) & 0x0F0F0F0F0F0F0F0FL;
      h = ((h <<  2) | h) & 0x3333333333333333L;
      h = ((h <<  1) | h) & 0x5555555555555555L;
      return h;
    }

    @Override
    public final long hash2ij(long h) {
      h = (0x2222222222222222L & h) <<  1
          | (0x4444444444444444L & h) >>> 1
          | (0x9999999999999999L & h);
      h = (0x0C0C0C0C0C0C0C0CL & h) <<  2
          | (0x3030303030303030L & h) >>> 2
          | (0xC3C3C3C3C3C3C3C3L & h);
      h = (0x00F000F000F000F0L & h) <<  4
          | (0x0F000F000F000F00L & h) >>> 4
          | (0xF00FF00FF00FF00FL & h);
      h = (0x0000FF000000FF00L & h) <<  8
          | (0x00FF000000FF0000L & h) >>> 8
          | (0xFF0000FFFF0000FFL & h);
      h = (0x00000000FFFF0000L & h) <<  16
          | (0x0000FFFF00000000L & h) >>> 16
          | (0xFFFF00000000FFFFL & h);
      return h;
    }

    @Override
    public final long hash2i0(long h) {
      // Verify all bit of j are set to 0
      assert (0x3333333333333333L & h) == 0;
      h = ((h >> 1) | h) & 0x3333333333333333L;
      h = ((h >> 2) | h) & 0x0F0F0F0F0F0F0F0FL;
      h = ((h >> 4) | h) & 0x00FF00FF00FF00FFL;
      h = ((h >> 8) | h) & 0x0000FFFF0000FFFFL;
      return h;
    }

    @Override
    public final int ij2i(final long ij) {
      return (int) ij;
    }

    @Override
    public final int ij2j(final long ij) {
      return (int) (ij >>> 32);
    }
  };

  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on the bitwise XOR operator to interleave the bits of
   * the discretized 2d-coordinates. We assume that each discritized coordinates is coded on maximum
   * 8 bits (BYTE).
   * The algorithm is a slightly adapted version of the outer perfect shuffle define
   * p. 106 of "Hacker's Delight" (Henry S. Warren, Jr), slightly modified.
   */
  public static final FillingCurve2D ZOC_VMSB_XOR_BYTE = new FillingCurve2D() {

    @Override
    public final long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public final long ij2hash(int i, int j) {
      // Verify for each coordinate that no bit is set after the 8th LSB
      assert (i== -1 && j== -1) || ((i & 0xFFFFFF00) == 0 && (j & 0xFFFFFF00) == 0);
      i |= (j << 8);
      j = (i ^ (i >>> 4)) & 0x000000F0; i = i ^ j ^ (j << 4);
      j = (i ^ (i >>> 2)) & 0x00000C0C; i = i ^ j ^ (j << 2);
      j = (i ^ (i >>> 1)) & 0x00002222; i = i ^ j ^ (j << 1);
      return i;
    }

    @Override
    public final long hash2ij(long hash) {
      // Verify no bit set after the 16th LSB
      assert (0xFFFFFFFFFFFF0000L & hash) == 0;
      int t, h = (int) hash;
      t = (h ^ (h >>> 1)) & 0x00002222; h = h ^ t ^ (t << 1);
      t = (h ^ (h >>> 2)) & 0x00000C0C; h = h ^ t ^ (t << 2);
      t = (h ^ (h >>> 4)) & 0x000000F0; h = h ^ t ^ (t << 4);
      return h;
    }

    @Override
    public final int ij2i(final long ij) {
      return (((int) ij) & 0x000000FF);
    }

    @Override
    public final int ij2j(final long ij) {
      return (((int) ij) >>> 8);
    }

    @Override
    public final long i02hash(final int i) {
      return ZOC_VMSB_OR_BYTE.i02hash(i);
    }

    @Override
    public final long hash2i0(final long hash) {
      return ZOC_VMSB_OR_BYTE.hash2i0(hash);
    }
  };

  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on the bitwise XOR operator to interleave the bits of
   * the discretized 2d-coordinates. We assume that each discritized coordinates is coded on maximum
   * 16 bits (SHORT).
   * The algorithm is a slightly adapted version of the outer perfect shuffle define
   * p. 106 of "Hacker's Delight" (Henry S. Warren, Jr).
   */
  public static final FillingCurve2D ZOC_VMSB_XOR_SHORT = new FillingCurve2D() {

    @Override
    public final long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public final long ij2hash(int i, int j) {
      // Verify for each coordinate that no bit is set after the 16th LSB
      assert (i == -1 && j == -1) || ((i & 0xFFFF0000) == 0 && (j & 0xFFFF0000) == 0);
      long k = (long) j;
      long h = ((long) i) | (k << 16);
      k = (h ^ (h >>> 8)) & 0x000000000000FF00L; h = h ^ k ^ (k << 8);
      k = (h ^ (h >>> 4)) & 0x0000000000F000F0L; h = h ^ k ^ (k << 4);
      k = (h ^ (h >>> 2)) & 0x000000000C0C0C0CL; h = h ^ k ^ (k << 2);
      k = (h ^ (h >>> 1)) & 0x0000000022222222L; h = h ^ k ^ (k << 1);
      return h;
    }

    @Override
    public final long hash2ij(long hash) {
      // Verify no bit set after the 32th LSB
      assert (0xFFFFFFFF00000000L & hash) == 0;
      int t, h = (int) hash;
      t = (h ^ (h >>> 1)) & 0x22222222; h = h ^ t ^ (t << 1);
      t = (h ^ (h >>> 2)) & 0x0C0C0C0C; h = h ^ t ^ (t << 2);
      t = (h ^ (h >>> 4)) & 0x00F000F0; h = h ^ t ^ (t << 4);
      t = (h ^ (h >>> 8)) & 0x0000FF00; h = h ^ t ^ (t << 8);
      return h;
    }

    @Override
    public final int ij2i(final long ij) {
      return (((int) ij) & 0x0000FFFF);
    }

    @Override
    public final int ij2j(final long ij) {
      return (((int) ij) >>> 16);
    }

    @Override
    public final long i02hash(final int i) {
      return ZOC_VMSB_OR_SHORT.i02hash(i);
    }

    @Override
    public final long hash2i0(final long hash) {
      return ZOC_VMSB_OR_SHORT.hash2i0(hash);
    }
  };


  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on the bitwise XOR operator to interleave the bits of
   * the discretized 2d-coordinates. We assume that each discritized coordinates is coded on maximum
   * 32 bits (INT).
   * The algorithm is a slightly adapted version of the outer perfect shuffle define
   * p. 106 of "Hacker's Delight" (Henry S. Warren, Jr), slightly modified.
   */
  public static final FillingCurve2D ZOC_VMSB_XOR_INT = new FillingCurve2D() {

    @Override
    public final long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public final long ij2hash(int i, int j) {
      long h = ((long) j) << 32, t;
      h |= i;
      t = (h ^ (h >> 16)) & 0x00000000FFFF0000L; h = h ^ t ^ (t << 16);
      t = (h ^ (h >>  8)) & 0x0000FF000000FF00L; h = h ^ t ^ (t <<  8);
      t = (h ^ (h >>  4)) & 0x00F000F000F000F0L; h = h ^ t ^ (t <<  4);
      t = (h ^ (h >>  2)) & 0x0C0C0C0C0C0C0C0CL; h = h ^ t ^ (t <<  2);
      t = (h ^ (h >>  1)) & 0x2222222222222222L; h = h ^ t ^ (t <<  1);
      return h;
    }

    @Override
    public final long i02hash(final int i) {
      return ZOC_VMSB_OR_INT.i02hash(i);
    }

    @Override
    public final long hash2ij(long h) {
      long t;
      t = (h ^ (h >>>  1)) & 0x2222222222222222L; h = h ^ t ^ (t <<  1);
      t = (h ^ (h >>>  2)) & 0x0C0C0C0C0C0C0C0CL; h = h ^ t ^ (t <<  2);
      t = (h ^ (h >>>  4)) & 0x00F000F000F000F0L; h = h ^ t ^ (t <<  4);
      t = (h ^ (h >>>  8)) & 0x0000FF000000FF00L; h = h ^ t ^ (t <<  8);
      t = (h ^ (h >>> 16)) & 0x00000000FFFF0000L; h = h ^ t ^ (t << 16);
      return h;
    }

    @Override
    public final long hash2i0(final long hash) {
      return ZOC_VMSB_OR_INT.hash2i0(hash);
    }

    @Override
    public final int ij2i(final long ij) {
      return (int) ij;
    }

    @Override
    public final int ij2j(final long ij) {
      return (int) (ij >>> 32);
    }
  };

  /** Lookup table storing the result of the {@link FillingCuve2D#i02hash(int)} method
   *  for all all possible unsigned bytes, so for all value in [0, 255].
   *  We tested the difference of performance between array of short and of long.
   *  In the first case, casts (into long) are required but the array is smaller so will occupy
   *  less space in the L1 cache. We were not able to detect any difference so we choose short.  */
  private static final short[] LUPT_TO_HASH = new short[] {
      0x0000, 0x0001, 0x0004, 0x0005, 0x0010, 0x0011, 0x0014, 0x0015, 0x0040, 0x0041, 0x0044,
      0x0045, 0x0050, 0x0051, 0x0054, 0x0055, 0x0100, 0x0101, 0x0104, 0x0105, 0x0110, 0x0111,
      0x0114, 0x0115, 0x0140, 0x0141, 0x0144, 0x0145, 0x0150, 0x0151, 0x0154, 0x0155, 0x0400,
      0x0401, 0x0404, 0x0405, 0x0410, 0x0411, 0x0414, 0x0415, 0x0440, 0x0441, 0x0444, 0x0445,
      0x0450, 0x0451, 0x0454, 0x0455, 0x0500, 0x0501, 0x0504, 0x0505, 0x0510, 0x0511, 0x0514,
      0x0515, 0x0540, 0x0541, 0x0544, 0x0545, 0x0550, 0x0551, 0x0554, 0x0555, 0x1000, 0x1001,
      0x1004, 0x1005, 0x1010, 0x1011, 0x1014, 0x1015, 0x1040, 0x1041, 0x1044, 0x1045, 0x1050,
      0x1051, 0x1054, 0x1055, 0x1100, 0x1101, 0x1104, 0x1105, 0x1110, 0x1111, 0x1114, 0x1115,
      0x1140, 0x1141, 0x1144, 0x1145, 0x1150, 0x1151, 0x1154, 0x1155, 0x1400, 0x1401, 0x1404,
      0x1405, 0x1410, 0x1411, 0x1414, 0x1415, 0x1440, 0x1441, 0x1444, 0x1445, 0x1450, 0x1451,
      0x1454, 0x1455, 0x1500, 0x1501, 0x1504, 0x1505, 0x1510, 0x1511, 0x1514, 0x1515, 0x1540,
      0x1541, 0x1544, 0x1545, 0x1550, 0x1551, 0x1554, 0x1555, 0x4000, 0x4001, 0x4004, 0x4005,
      0x4010, 0x4011, 0x4014, 0x4015, 0x4040, 0x4041, 0x4044, 0x4045, 0x4050, 0x4051, 0x4054,
      0x4055, 0x4100, 0x4101, 0x4104, 0x4105, 0x4110, 0x4111, 0x4114, 0x4115, 0x4140, 0x4141,
      0x4144, 0x4145, 0x4150, 0x4151, 0x4154, 0x4155, 0x4400, 0x4401, 0x4404, 0x4405, 0x4410,
      0x4411, 0x4414, 0x4415, 0x4440, 0x4441, 0x4444, 0x4445, 0x4450, 0x4451, 0x4454, 0x4455,
      0x4500, 0x4501, 0x4504, 0x4505, 0x4510, 0x4511, 0x4514, 0x4515, 0x4540, 0x4541, 0x4544,
      0x4545, 0x4550, 0x4551, 0x4554, 0x4555, 0x5000, 0x5001, 0x5004, 0x5005, 0x5010, 0x5011,
      0x5014, 0x5015, 0x5040, 0x5041, 0x5044, 0x5045, 0x5050, 0x5051, 0x5054, 0x5055, 0x5100,
      0x5101, 0x5104, 0x5105, 0x5110, 0x5111, 0x5114, 0x5115, 0x5140, 0x5141, 0x5144, 0x5145,
      0x5150, 0x5151, 0x5154, 0x5155, 0x5400, 0x5401, 0x5404, 0x5405, 0x5410, 0x5411, 0x5414,
      0x5415, 0x5440, 0x5441, 0x5444, 0x5445, 0x5450, 0x5451, 0x5454, 0x5455, 0x5500, 0x5501,
      0x5504, 0x5505, 0x5510, 0x5511, 0x5514, 0x5515, 0x5540, 0x5541, 0x5544, 0x5545, 0x5550,
      0x5551, 0x5554, 0x5555
  };

  private static final short[] LUPT_TO_IJ_BYTE = new short[] {
      0x000, 0x001, 0x100, 0x101, 0x002, 0x003, 0x102, 0x103, 0x200, 0x201, 0x300, 0x301, 0x202,
      0x203, 0x302, 0x303, 0x004, 0x005, 0x104, 0x105, 0x006, 0x007, 0x106, 0x107, 0x204, 0x205,
      0x304, 0x305, 0x206, 0x207, 0x306, 0x307, 0x400, 0x401, 0x500, 0x501, 0x402, 0x403, 0x502,
      0x503, 0x600, 0x601, 0x700, 0x701, 0x602, 0x603, 0x702, 0x703, 0x404, 0x405, 0x504, 0x505,
      0x406, 0x407, 0x506, 0x507, 0x604, 0x605, 0x704, 0x705, 0x606, 0x607, 0x706, 0x707, 0x008,
      0x009, 0x108, 0x109, 0x00A, 0x00B, 0x10A, 0x10B, 0x208, 0x209, 0x308, 0x309, 0x20A, 0x20B,
      0x30A, 0x30B, 0x00C, 0x00D, 0x10C, 0x10D, 0x00E, 0x00F, 0x10E, 0x10F, 0x20C, 0x20D, 0x30C,
      0x30D, 0x20E, 0x20F, 0x30E, 0x30F, 0x408, 0x409, 0x508, 0x509, 0x40A, 0x40B, 0x50A, 0x50B,
      0x608, 0x609, 0x708, 0x709, 0x60A, 0x60B, 0x70A, 0x70B, 0x40C, 0x40D, 0x50C, 0x50D, 0x40E,
      0x40F, 0x50E, 0x50F, 0x60C, 0x60D, 0x70C, 0x70D, 0x60E, 0x60F, 0x70E, 0x70F, 0x800, 0x801,
      0x900, 0x901, 0x802, 0x803, 0x902, 0x903, 0xA00, 0xA01, 0xB00, 0xB01, 0xA02, 0xA03, 0xB02,
      0xB03, 0x804, 0x805, 0x904, 0x905, 0x806, 0x807, 0x906, 0x907, 0xA04, 0xA05, 0xB04, 0xB05,
      0xA06, 0xA07, 0xB06, 0xB07, 0xC00, 0xC01, 0xD00, 0xD01, 0xC02, 0xC03, 0xD02, 0xD03, 0xE00,
      0xE01, 0xF00, 0xF01, 0xE02, 0xE03, 0xF02, 0xF03, 0xC04, 0xC05, 0xD04, 0xD05, 0xC06, 0xC07,
      0xD06, 0xD07, 0xE04, 0xE05, 0xF04, 0xF05, 0xE06, 0xE07, 0xF06, 0xF07, 0x808, 0x809, 0x908,
      0x909, 0x80A, 0x80B, 0x90A, 0x90B, 0xA08, 0xA09, 0xB08, 0xB09, 0xA0A, 0xA0B, 0xB0A, 0xB0B,
      0x80C, 0x80D, 0x90C, 0x90D, 0x80E, 0x80F, 0x90E, 0x90F, 0xA0C, 0xA0D, 0xB0C, 0xB0D, 0xA0E,
      0xA0F, 0xB0E, 0xB0F, 0xC08, 0xC09, 0xD08, 0xD09, 0xC0A, 0xC0B, 0xD0A, 0xD0B, 0xE08, 0xE09,
      0xF08, 0xF09, 0xE0A, 0xE0B, 0xF0A, 0xF0B, 0xC0C, 0xC0D, 0xD0C, 0xD0D, 0xC0E, 0xC0F, 0xD0E,
      0xD0F, 0xE0C, 0xE0D, 0xF0C, 0xF0D, 0xE0E, 0xE0F, 0xF0E, 0xF0F
  };

  private static final int[] LUPT_TO_IJ_SHORT = new int[] {
      0x00000, 0x00001, 0x10000, 0x10001, 0x00002, 0x00003, 0x10002, 0x10003, 0x20000, 0x20001,
      0x30000, 0x30001, 0x20002, 0x20003, 0x30002, 0x30003, 0x00004, 0x00005, 0x10004, 0x10005,
      0x00006, 0x00007, 0x10006, 0x10007, 0x20004, 0x20005, 0x30004, 0x30005, 0x20006, 0x20007,
      0x30006, 0x30007, 0x40000, 0x40001, 0x50000, 0x50001, 0x40002, 0x40003, 0x50002, 0x50003,
      0x60000, 0x60001, 0x70000, 0x70001, 0x60002, 0x60003, 0x70002, 0x70003, 0x40004, 0x40005,
      0x50004, 0x50005, 0x40006, 0x40007, 0x50006, 0x50007, 0x60004, 0x60005, 0x70004, 0x70005,
      0x60006, 0x60007, 0x70006, 0x70007, 0x00008, 0x00009, 0x10008, 0x10009, 0x0000A, 0x0000B,
      0x1000A, 0x1000B, 0x20008, 0x20009, 0x30008, 0x30009, 0x2000A, 0x2000B, 0x3000A, 0x3000B,
      0x0000C, 0x0000D, 0x1000C, 0x1000D, 0x0000E, 0x0000F, 0x1000E, 0x1000F, 0x2000C, 0x2000D,
      0x3000C, 0x3000D, 0x2000E, 0x2000F, 0x3000E, 0x3000F, 0x40008, 0x40009, 0x50008, 0x50009,
      0x4000A, 0x4000B, 0x5000A, 0x5000B, 0x60008, 0x60009, 0x70008, 0x70009, 0x6000A, 0x6000B,
      0x7000A, 0x7000B, 0x4000C, 0x4000D, 0x5000C, 0x5000D, 0x4000E, 0x4000F, 0x5000E, 0x5000F,
      0x6000C, 0x6000D, 0x7000C, 0x7000D, 0x6000E, 0x6000F, 0x7000E, 0x7000F, 0x80000, 0x80001,
      0x90000, 0x90001, 0x80002, 0x80003, 0x90002, 0x90003, 0xA0000, 0xA0001, 0xB0000, 0xB0001,
      0xA0002, 0xA0003, 0xB0002, 0xB0003, 0x80004, 0x80005, 0x90004, 0x90005, 0x80006, 0x80007,
      0x90006, 0x90007, 0xA0004, 0xA0005, 0xB0004, 0xB0005, 0xA0006, 0xA0007, 0xB0006, 0xB0007,
      0xC0000, 0xC0001, 0xD0000, 0xD0001, 0xC0002, 0xC0003, 0xD0002, 0xD0003, 0xE0000, 0xE0001,
      0xF0000, 0xF0001, 0xE0002, 0xE0003, 0xF0002, 0xF0003, 0xC0004, 0xC0005, 0xD0004, 0xD0005,
      0xC0006, 0xC0007, 0xD0006, 0xD0007, 0xE0004, 0xE0005, 0xF0004, 0xF0005, 0xE0006, 0xE0007,
      0xF0006, 0xF0007, 0x80008, 0x80009, 0x90008, 0x90009, 0x8000A, 0x8000B, 0x9000A, 0x9000B,
      0xA0008, 0xA0009, 0xB0008, 0xB0009, 0xA000A, 0xA000B, 0xB000A, 0xB000B, 0x8000C, 0x8000D,
      0x9000C, 0x9000D, 0x8000E, 0x8000F, 0x9000E, 0x9000F, 0xA000C, 0xA000D, 0xB000C, 0xB000D,
      0xA000E, 0xA000F, 0xB000E, 0xB000F, 0xC0008, 0xC0009, 0xD0008, 0xD0009, 0xC000A, 0xC000B,
      0xD000A, 0xD000B, 0xE0008, 0xE0009, 0xF0008, 0xF0009, 0xE000A, 0xE000B, 0xF000A, 0xF000B,
      0xC000C, 0xC000D, 0xD000C, 0xD000D, 0xC000E, 0xC000F, 0xD000E, 0xD000F, 0xE000C, 0xE000D,
      0xF000C, 0xF000D, 0xE000E, 0xE000F, 0xF000E, 0xF000F
  };

  private static final long[] LUPT_TO_IJ_INT = new long[] {
      0x000000000L, 0x000000001L, 0x100000000L, 0x100000001L, 0x000000002L, 0x000000003L,
      0x100000002L, 0x100000003L, 0x200000000L, 0x200000001L, 0x300000000L, 0x300000001L,
      0x200000002L, 0x200000003L, 0x300000002L, 0x300000003L, 0x000000004L, 0x000000005L,
      0x100000004L, 0x100000005L, 0x000000006L, 0x000000007L, 0x100000006L, 0x100000007L,
      0x200000004L, 0x200000005L, 0x300000004L, 0x300000005L, 0x200000006L, 0x200000007L,
      0x300000006L, 0x300000007L, 0x400000000L, 0x400000001L, 0x500000000L, 0x500000001L,
      0x400000002L, 0x400000003L, 0x500000002L, 0x500000003L, 0x600000000L, 0x600000001L,
      0x700000000L, 0x700000001L, 0x600000002L, 0x600000003L, 0x700000002L, 0x700000003L,
      0x400000004L, 0x400000005L, 0x500000004L, 0x500000005L, 0x400000006L, 0x400000007L,
      0x500000006L, 0x500000007L, 0x600000004L, 0x600000005L, 0x700000004L, 0x700000005L,
      0x600000006L, 0x600000007L, 0x700000006L, 0x700000007L, 0x000000008L, 0x000000009L,
      0x100000008L, 0x100000009L, 0x00000000AL, 0x00000000BL, 0x10000000AL, 0x10000000BL,
      0x200000008L, 0x200000009L, 0x300000008L, 0x300000009L, 0x20000000AL, 0x20000000BL,
      0x30000000AL, 0x30000000BL, 0x00000000CL, 0x00000000DL, 0x10000000CL, 0x10000000DL,
      0x00000000EL, 0x00000000FL, 0x10000000EL, 0x10000000FL, 0x20000000CL, 0x20000000DL,
      0x30000000CL, 0x30000000DL, 0x20000000EL, 0x20000000FL, 0x30000000EL, 0x30000000FL,
      0x400000008L, 0x400000009L, 0x500000008L, 0x500000009L, 0x40000000AL, 0x40000000BL,
      0x50000000AL, 0x50000000BL, 0x600000008L, 0x600000009L, 0x700000008L, 0x700000009L,
      0x60000000AL, 0x60000000BL, 0x70000000AL, 0x70000000BL, 0x40000000CL, 0x40000000DL,
      0x50000000CL, 0x50000000DL, 0x40000000EL, 0x40000000FL, 0x50000000EL, 0x50000000FL,
      0x60000000CL, 0x60000000DL, 0x70000000CL, 0x70000000DL, 0x60000000EL, 0x60000000FL,
      0x70000000EL, 0x70000000FL, 0x800000000L, 0x800000001L, 0x900000000L, 0x900000001L,
      0x800000002L, 0x800000003L, 0x900000002L, 0x900000003L, 0xA00000000L, 0xA00000001L,
      0xB00000000L, 0xB00000001L, 0xA00000002L, 0xA00000003L, 0xB00000002L, 0xB00000003L,
      0x800000004L, 0x800000005L, 0x900000004L, 0x900000005L, 0x800000006L, 0x800000007L,
      0x900000006L, 0x900000007L, 0xA00000004L, 0xA00000005L, 0xB00000004L, 0xB00000005L,
      0xA00000006L, 0xA00000007L, 0xB00000006L, 0xB00000007L, 0xC00000000L, 0xC00000001L,
      0xD00000000L, 0xD00000001L, 0xC00000002L, 0xC00000003L, 0xD00000002L, 0xD00000003L,
      0xE00000000L, 0xE00000001L, 0xF00000000L, 0xF00000001L, 0xE00000002L, 0xE00000003L,
      0xF00000002L, 0xF00000003L, 0xC00000004L, 0xC00000005L, 0xD00000004L, 0xD00000005L,
      0xC00000006L, 0xC00000007L, 0xD00000006L, 0xD00000007L, 0xE00000004L, 0xE00000005L,
      0xF00000004L, 0xF00000005L, 0xE00000006L, 0xE00000007L, 0xF00000006L, 0xF00000007L,
      0x800000008L, 0x800000009L, 0x900000008L, 0x900000009L, 0x80000000AL, 0x80000000BL,
      0x90000000AL, 0x90000000BL, 0xA00000008L, 0xA00000009L, 0xB00000008L, 0xB00000009L,
      0xA0000000AL, 0xA0000000BL, 0xB0000000AL, 0xB0000000BL, 0x80000000CL, 0x80000000DL,
      0x90000000CL, 0x90000000DL, 0x80000000EL, 0x80000000FL, 0x90000000EL, 0x90000000FL,
      0xA0000000CL, 0xA0000000DL, 0xB0000000CL, 0xB0000000DL, 0xA0000000EL, 0xA0000000FL,
      0xB0000000EL, 0xB0000000FL, 0xC00000008L, 0xC00000009L, 0xD00000008L, 0xD00000009L,
      0xC0000000AL, 0xC0000000BL, 0xD0000000AL, 0xD0000000BL, 0xE00000008L, 0xE00000009L,
      0xF00000008L, 0xF00000009L, 0xE0000000AL, 0xE0000000BL, 0xF0000000AL, 0xF0000000BL,
      0xC0000000CL, 0xC0000000DL, 0xD0000000CL, 0xD0000000DL, 0xC0000000EL, 0xC0000000FL,
      0xD0000000EL, 0xD0000000FL, 0xE0000000CL, 0xE0000000DL, 0xF0000000CL, 0xF0000000DL,
      0xE0000000EL, 0xE0000000FL, 0xF0000000EL, 0xF0000000FL
  };

  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on a lookup table (LOOKUP).
   * We assume that each discritized coordinates is coded on maximum 8 bits (BYTE).
   */
  public static final FillingCurve2D ZOC_VMSB_LOOKUP_BYTE = new FillingCurve2D() {

    @Override
    public long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public long ij2hash(int i, int j) {
      return (i02hash(j) << 1) | i02hash(i);
    }

    @Override
    public long i02hash(int i) {
      // Verify that the coordinate has no bit set after the 8th LSB
      assert (i == -1) || (i & 0xFFFFFF00) == 0;
      return LUPT_TO_HASH[i & 0x00000FF]; // apply the mask just in case of -1 value
    }

    @Override
    public long hash2ij(long hash) {
      // Verify no bit set after the 16th LSB
      assert (0xFFFFFFFFFFFF0000L & hash) == 0;
      int h = (int) hash;
      return LUPT_TO_IJ_BYTE[h >> 8] << 4
          | LUPT_TO_IJ_BYTE[h & 0x000000FF];
    }

    @Override
    public long hash2i0(long hash) {
      // Verify all bit of j are set to 0
      assert (0xFFFFFFFFFFFFAAAAL & hash) == 0;
      return hash2ij(hash);
    }

    @Override
    public int ij2i(final long ij) {
      return (((int) ij) & 0x000000FF);
    }

    @Override
    public int ij2j(final long ij) {
      return (((int) ij) >>> 8);
    }
  };


  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on a lookup table (LOOKUP).
   * We assume that each discritized coordinates is coded on maximum 16 bits (SHORT).
   */
  public static final FillingCurve2D ZOC_VMSB_LOOKUP_SHORT = new FillingCurve2D() {

    @Override
    public long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public long ij2hash(int i, int j) {
       return (i02hash(j) << 1) | i02hash(i);
    }

    @Override
    public long i02hash(int i) {
      // Verify that the coordinate has no bit set after the 8th LSB
      assert (i == -1) || ((i & 0xFFFF0000) == 0);
      return ((long) LUPT_TO_HASH[(i >>> 8) & 0x000000FF]) << 16 // apply the mask just in case of -1 value
          | LUPT_TO_HASH[i & 0x000000FF];
    }

    @Override
    public long hash2ij(long hash) {
      // Verify no bit set after the 16th LSB
      assert (0xFFFFFFFF00000000L & hash) == 0;
      return LUPT_TO_IJ_SHORT[(int) ((hash & 0x00000000FF000000L) >> 24) ] << 12
           | LUPT_TO_IJ_SHORT[(int) ((hash & 0x0000000000FF0000L) >> 16) ] <<  8
           | LUPT_TO_IJ_SHORT[(int) ((hash & 0x000000000000FF00L) >>  8) ] <<  4
           | LUPT_TO_IJ_SHORT[(int)  (hash & 0x00000000000000FFL)];

    }

    @Override
    public long hash2i0(long hash) {
      // Verify all bit of j are set to 0
      assert (0xFFFFFFFFAAAAAAAAL & hash) == 0;
      return hash2ij(hash);
    }

    @Override
    public int ij2i(final long ij) {
      return (((int) ij) & 0x0000FFFF);
    }

    @Override
    public int ij2j(final long ij) {
      return (((int) ij) >>> 16);
    }
  };

  /**
   * Z-Order Curve (ZOC) implementation in which the vertical coordinate carry the most significant
   * bit (VMSB). This implementation is based on a lookup table (LOOKUP).
   * We assume that each discritized coordinates is coded on maximum 32 bits (INT).
   */
  public static final FillingCurve2D ZOC_VMSB_LOOKUP_INT = new FillingCurve2D() {

    @Override
    public long xy2hash(final double x, final double y) {
      return ij2hash((int) x, (int) y);
    }

    @Override
    public long ij2hash(int i, int j) {
      return (i02hash(j) << 1) | i02hash(i);
    }

    @Override
    public long i02hash(int i) {
      // (long) (LUPT_TO_HASH[(i & 0xFF000000) >>> 24]) << 48
      return (long) (LUPT_TO_HASH[i >>> 24]) << 48
           | (long) (LUPT_TO_HASH[(i & 0x00FF0000) >>> 16]) << 32
           | (long) (LUPT_TO_HASH[(i & 0x0000FF00) >>>  8]) << 16
           | (long) (LUPT_TO_HASH[ i & 0x000000FF]);
    }

    @Override
    public long hash2ij(long h) {
      return LUPT_TO_IJ_INT[(int) ((h & 0xFF00000000000000L) >>> 56)] << 28
           | LUPT_TO_IJ_INT[(int) ((h & 0x00FF000000000000L) >>> 48)] << 24
           | LUPT_TO_IJ_INT[(int) ((h & 0x0000FF0000000000L) >>> 40)] << 20
           | LUPT_TO_IJ_INT[(int) ((h & 0x000000FF00000000L) >>> 32)] << 16
           | LUPT_TO_IJ_INT[(int) ((h & 0x00000000FF000000L) >>> 24)] << 12
           | LUPT_TO_IJ_INT[(int) ((h & 0x0000000000FF0000L) >>> 16)] <<  8
           | LUPT_TO_IJ_INT[(int) ((h & 0x000000000000FF00L) >>>  8)] <<  4
           | LUPT_TO_IJ_INT[(int)  (h & 0x00000000000000FFL)];
    }

    @Override
    public long hash2i0(long hash) {
      // Verify all bit of j are set to 0
      assert (0xFFFFFFFF33333333L & hash) == 0;
      return hash2ij(hash);
    }

    @Override
    public final int ij2i(final long ij) {
      return (int) ij;
    }

    @Override
    public final int ij2j(final long ij) {
      return (int) (ij >>> 32);
    }
  };

  /*private static final void buildLookUpTables() {
    for (int i = 0; i < 256; i++) {
      System.out.print(String.format("%0#6X", ZOC_VMSB_OR_BYTE.i02hash(i)) + ", ");
    }
    System.out.println("\n---- byte");
    for (int i = 0; i < 256; i++) {
      System.out.print(String.format("0x%03X", ZOC_VMSB_OR_BYTE.hash2ij(i)) + ", ");
    }
    System.out.println("\n---- short");
    for (int i = 0; i < 256; i++) {
      System.out.print(String.format("0x%05X", ZOC_VMSB_OR_SHORT.hash2ij(i)) + ", ");
    }
    System.out.println("\n---- int");
    for (int i = 0; i < 256; i++) {
      System.out.print(String.format("0x%09X", ZOC_VMSB_OR_INT.hash2ij(i)) + "L, ");
    }
  }
  public static void main(final String[] args) {
    buildLookUpTables();
  }*/


}
