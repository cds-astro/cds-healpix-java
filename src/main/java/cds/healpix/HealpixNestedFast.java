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

import java.util.Arrays;
import java.util.EnumMap;
import java.util.EnumSet;

import cds.healpix.fillingcurve.FillingCurve2D;
import cds.healpix.fillingcurve.FillingCurve2DType;
import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.CompassPoint.MainWind;
import cds.healpix.CompassPoint.Ordinal;
import cds.healpix.common.math.FastMath;

import static cds.healpix.Healpix.TRANSITION_LATITUDE;
import static cds.healpix.Healpix.TRANSITION_Z;
import static cds.healpix.Healpix.nside;
import static cds.healpix.common.math.FastMath.acos;
import static cds.healpix.common.math.FastMath.asin;
import static cds.healpix.common.math.HackersDelight.BUT_SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_I;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.toBits;
import static cds.healpix.common.math.Math.FOUR_OVER_PI;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.SQRT6;
import static cds.healpix.common.math.Math.pow;
import static cds.healpix.Healpix.ONE_OVER_TRANSITION_Z;
import static cds.healpix.Healpix.checkLatitude;
import static cds.healpix.HealpixUnprojector.ONE_OVER_SQRT6;
import static cds.healpix.HealpixNested.bits2hash;

/**
 * A faster, thread-safe, but ugly (sorry: less readable) version of {@link HealpixNested}.
 * 
 * The code has also been "denormalise" to try to be faster: we introduce redundancy by suppressing 
 * method calls when a method must return more than one elements, ...
 * 
 * Motivations: in Java, all arguments of a function are passed-by-value. And, like most languages,
 * a function returns only one output. It means that if, e.g. a function as to return a pair of
 * coordinate you have only two solutions: the first one is to instantiate an object containing
 * the pair of coordinates and return it, or the method have to take an object and set the pair 
 * of coordinates.
 * In this "fast version", we duplicate code to avoid creating objects inside critical methods and
 * we write the results in arguments instead of returning new objects. The instantiation of an
 * object in Java is cheap, but it cost is not negligible compared to the speed of HEALPix methods. 
 * 
 * @author F.-X. Pineau
 *
 */
public final class HealpixNestedFast implements HashComputer, VerticesAndPathComputer, NeighbourSelector {

  /**
   * For each of the 12 hash of depth 0, contains the list of the depth 0
   * neighbours sorted by increasing hash value.
   */
  public static final byte[][] D0H_NEIGHBOURS = new byte[][]{
    // North polar cap (no East, no West)
    {1, 2, 3, 4, 5, 8},  // dh0 = 0
    {0, 2, 3, 5, 6, 9},  // dh0 = 1
    {0, 1, 3, 6, 7, 10}, // dh0 = 2
    {0, 1, 2, 4, 7, 11}, // dh0 = 3
    // Equatorial region (no North, not South)
    {0, 3, 5, 7, 8, 11}, // dh0 = 4
    {0, 1, 4, 6, 8, 9},  // dh0 = 5
    {1, 2, 5, 7, 9, 10}, // dh0 = 6
    {2, 3, 4, 6, 10, 11},// dh0 = 7
    // South polar cap (no East, no West)
    {0, 4, 5, 9, 10, 11},// dh0 = 8
    {1, 5, 6, 8, 10, 11},// dh0 = 9
    {2, 6, 7, 8, 9, 11}, // dh0 = 10
    {3, 4, 7, 8, 9, 10}  // dh0 = 11
  };

  public static final byte[][][] NEIGHBOURS = new byte[12][][];

  private final int depth;

  // Derived quantities to speed up computations
  private final int twiceDepth;
  private final int nside;
  private final int twiceNside;
  private final long halfNside4IEEEdouble;
  private final double halfNside;
  private final double oneOverNside;
  private final double piOverFourNside;
  private final int nsideRemainderMask;

  private final long nHash;

  /** Mask used to retrieve the "depth" least significant bits, i.e. this mask
   * is used to compute "i modulo 2^depth", i.e. "i modulo nside". */
  private final int modNsideMask;
  /** Mask used to retrieve the value of the base (i.e., depth=0) pixel in
   * the full hash. So all bits are : 15 << (this.depth << 1).
   * Background: we code 2^n - 1 values on n bits */
  private final long d0Mask;
  /** Mask used to retrieve the yx part of a hash, i.e. removing the value
   * of the base pixel: the 2*nside least significant bits are set to 1.  */
  private final long xyMask;
  /** Mask used to retrieve the bits of the x value only in the full hash. */
  private final long xMask;
  /** Mask used to retrieve the bits of the y value only in the full hash. */
  private final long yMask;

  private final FillingCurve2D fc;

  public HealpixNestedFast(final int depth, FillingCurve2DType fillingCurveType) {
    this.depth = depth;
    this.nside = nside(this.depth);
    this.twiceNside = this.nside << 1;
    // this.nside4IEEEdouble = Healpix.nside4IEEEdouble(depth);
    this.halfNside4IEEEdouble = Healpix.halfNside4IEEEdouble(depth);
    this.halfNside = this.nside >> 1;
    this.oneOverNside = 1 / (double) nside;
    this.piOverFourNside = PI_OVER_FOUR * this.oneOverNside;
    this.twiceDepth = this.depth << 1;
    this.nsideRemainderMask = this.nside - 1;

    this.nHash = Healpix.nHash(this.depth);

    this.d0Mask = 15L << this.twiceDepth; // 15: 0...01111 in binary
    this.xyMask = (1L << this.twiceDepth) - 1;
    this.xMask = 0x5555555555555555L >>> (64 - this.twiceDepth); // ...0101
    this.yMask = 0xAAAAAAAAAAAAAAAAL >>> (64 - this.twiceDepth); // ...1010 = xMask << 1
    this.modNsideMask = this.nside - 1;
    // Set the filling curve object to be used
    this.fc = fillingCurveType.get(depth);
  }

  @Override
  public int depth() {
    return this.depth;
  }

  /////////////////////////////////
  // HashComputer implementation //
  /////////////////////////////////

  @Override
  public long hash(double lonRad, double latRad) {
    // We change the algo with respect to the cleaner code in HealpixNestedHashComputer.
    // The reasons are:
    // - see class Javadoc (for performacnes: we introduce redundancy, inline all method calls, ...)
    // - as we introduce redondancy, we do not have to keep the projection code separated
    //   - it allows to perform the computations in the 2 canonical base cells (NPC: 0; EQR: 4),
    //     with the y-origin of the Collignon projection set to 0 (instead of 1), and only then 
    //     we consider the base cell number
    // - thus we also removed the base cell lookup table:
    //   - it saves CPU cache utilisation and lookup
    //   - but it introduces 'if' conditions that we tried to avoid: effect should be limited since 
    //     few risks of branche miss-predictions, except in South/North polar caps (effect should be
    //     negligeable, except if input positions are random and more or less uniformly distributed
    //     in South and North polar caps).
    checkLatitude(latRad);
    final double absLon = fromBits(toBits(lonRad) & BUT_SIGN_BIT_MASK_L);  assert  0 <= absLon && absLon <= 10 * TWO_PI : absLon; // 10 is arbitrary
    long signLat = toBits(latRad);                                         assert latRad == fromBits(signLat);
    final double absLat = fromBits(signLat & BUT_SIGN_BIT_MASK_L);         assert  0 <= absLat && absLat <= HALF_PI; 
    signLat &= SIGN_BIT_MASK_L;                                            assert signLat == 0 || signLat == SIGN_BIT_MASK_L;
    double x = FOUR_OVER_PI * absLon, y;
    int xfloor = ((int) x);
    x -= (xfloor | 1);                                                     assert -1 <= x && x <= 1;
    long d0h;
    int i, j;
    if (absLat <= TRANSITION_LATITUDE) { // Equatorial region (Cylindrical projection)
      y = FastMath.sinQ(absLat, 0);                                        assert  0 <= y && y <= TRANSITION_Z;
      y *= ONE_OVER_TRANSITION_Z;                                          assert  0 <= 0 && y <= 1;
      y = fromBits(signLat | toBits(y));                                   assert -1 <= y && y <= 1;
      y += 2;                                                              assert  1 <= y && y <= 3;
      // Rotation of 45 and scale by sqrt(2) * nside / 2
      i = (int) fromBits(this.halfNside4IEEEdouble + toBits(y + x));       assert 0 <= i && i <= 2 * nside;
      if (i == twiceNside) { --i; }    // if, but very rare case of branch miss-prediction
      j = (int) fromBits(this.halfNside4IEEEdouble + toBits(y - x));       assert 0 <= j && j <= 2 * nside;
      if (j == twiceNside) { --j; }    // if, but very rare case of branch miss-prediction
      // Base cell hash computation
      final int id0c = i >>> this.depth;                                   assert id0c == 0 || id0c == 1;
      final int jd0c = j >>> this.depth;                                   assert jd0c == 0 || jd0c == 1;
      final int zeroOr1 = id0c ^ jd0c;                                     assert zeroOr1 == 0 || zeroOr1 == 1;
      final int zeroOr4Or8 = (2 >>> (id0c + jd0c)) << 2;                   assert zeroOr4Or8 == 0 || zeroOr4Or8 == 4 || zeroOr4Or8 == 8;
      d0h = zeroOr4Or8 | (((xfloor + zeroOr1) & 7) >> 1);                  assert 0 <= d0h && d0h < 12;
      // Coordinates in base cell
      i &= this.nsideRemainderMask;                                        assert 0 <= i && i < nside;
      j &= this.nsideRemainderMask;                                        assert 0 <= j && j < nside;
    } else { // Polar caps (Collignon projection)
      y = SQRT6 * FastMath.cosQ(0.5 * absLat + PI_OVER_FOUR);              assert  0 <= y && y <  1; 
      x *= y;                                                              assert -1 <= x && x <= 1;
      d0h = (xfloor & 7) >> 1;
      // Rotation of 45 and scale by sqrt(2) * nside / 2
      i = ((int) fromBits(this.halfNside4IEEEdouble + toBits(x + y)));     assert 0 <= i && i < nside;
      j = ((int) fromBits(this.halfNside4IEEEdouble + toBits(y - x)));     assert 0 <= j && j < nside;
      if (signLat == 0L) { // North polar cap
        xfloor = i;
        i = this.nsideRemainderMask - j;                                   assert 0 <= i && i < nside;
        j = this.nsideRemainderMask - xfloor;                              assert 0 <= j && j < nside;
      } else { // South polar cap
        d0h |= 8;                                                          assert 8 <= d0h && d0h < 12;
      }
    }
    return (d0h << this.twiceDepth) | this.fc.ij2hash(i, j);
  }

  ////////////////////////////////////////////
  // VerticesAndPathComputer implementation //
  ////////////////////////////////////////////

  @Override
  public double[] center(long hash) {
    final double[] resultLonLat = new double[2];
    center(hash, resultLonLat);
    return resultLonLat;
  }

  @Override
  public void center(long hash, final double[] resultLonLat) {
    checkHashRange(hash);
    // Pull apart the hash elements
    final int d0h = (int) (hash >> this.twiceDepth);         assert 0 <= d0h && d0h < 12;
    hash = this.fc.hash2ij(hash & this.xyMask);
    final int iInD0h = this.fc.ij2i(hash);                   assert 0 <= iInD0h && iInD0h < this.nside;
    final int jInD0h = this.fc.ij2j(hash);                   assert 0 <= jInD0h && jInD0h < this.nside;
    // Compute coordinates from the center of the base pixel with x-axis = W-->E, y-axis = S-->N
    final int lInD0h = iInD0h - jInD0h;                      assert -nside < lInD0h && lInD0h < nside;
    final int hInD0h = iInD0h + jInD0h - this.modNsideMask;  assert -nside < hInD0h && hInD0h < nside;
    // Compute coordinates of the base pixel in the projection plane
    final int d0hBy4Quotient = d0h >> 2;                     assert 0 <= d0hBy4Quotient && d0hBy4Quotient <= 2;
    final int d0hMod4 = d0h & 3;                             assert 0 <= d0hMod4 && d0hMod4 <= 3;
    final int hD0h = 1 - d0hBy4Quotient;                     assert -1 <= hD0h && hD0h <= 1;
    // +1 if the base cell is not equatorial
    int lD0h = (d0hMod4 << 1) | (hD0h & 1);                  assert (hD0h == 0 && lD0h == 0 || lD0h == 2 || lD0h == 4 ||  lD0h == 6)
                                                                 || (hD0h != 0 && lD0h == 1 || lD0h == 3 || lD0h == 5 ||  lD0h == 7);
    // Compute let's go
    double lon = lInD0h * this.oneOverNside;                  assert -1 < lon && lon < 1;
    double lat = hInD0h * this.oneOverNside;                  assert -1 < lat && lat < 1;
    if ((d0h < 4 && hInD0h > 0) || (d0h > 7 && hInD0h < 0)) { // Polar cap
      long signLat = toBits(lat);                            assert lat == fromBits(signLat);
      lat = fromBits(signLat & BUT_SIGN_BIT_MASK_L);         assert  0 < lat && lat < 1; 
      signLat &= SIGN_BIT_MASK_L;                            assert signLat == 0 || signLat == SIGN_BIT_MASK_L;
      lat = 1 - lat;                                         assert   0 < lat && lat <  1;
      lon /= lat;                                            assert  -1 < lon && lon <  1;
      lat *= ONE_OVER_SQRT6;                                 assert   0 < lat && lat < ONE_OVER_SQRT6;
      lat = 2 * acos(lat) - HALF_PI;                         assert TRANSITION_LATITUDE < lat && lat <= HALF_PI;
      lat = fromBits(signLat | toBits(lat));
    } else { // Equatorial region
      lat += hD0h;
      lat = asin(lat * TRANSITION_Z);
      // lD0h += 8 if d0h == 4 && lon < 0
      lD0h |= ((lInD0h & SIGN_BIT_MASK_I) >>> (24 + d0h)) & 8;  assert 0 <= lD0h && lD0h <= 8;
    }
    lon += lD0h;
    lon *= PI_OVER_FOUR;
    resultLonLat[LON_INDEX] = lon;
    resultLonLat[LAT_INDEX] = lat;
  }

  @Override
  public double[] vertex(long hash, Cardinal cardinalPoint) {
    final double[] resultLonLat = new double[2];
    vertex(hash, cardinalPoint, resultLonLat);
    return resultLonLat;
  }

  @Override
  public void vertex(long hash, Cardinal vertexDirection, double[] resultLonLat) {
    checkHashRange(hash);
    // Pull apart the hash elements
    final int d0h = (int) (hash >> this.twiceDepth);         assert 0 <= d0h && d0h < 12;
    hash = this.fc.hash2ij(hash & this.xyMask);
    final int iInD0h = this.fc.ij2i(hash);                   assert 0 <= iInD0h && iInD0h < this.nside;
    final int jInD0h = this.fc.ij2j(hash);                   assert 0 <= jInD0h && jInD0h < this.nside;
    // Compute coordinates from the center of the base pixel with x-axis = W-->E, y-axis = S-->N
    final int lInD0h = iInD0h - jInD0h;                      assert -nside < lInD0h && lInD0h < nside;
    final int hInD0h = iInD0h + jInD0h - this.modNsideMask;  assert -nside < hInD0h && hInD0h < nside;
    // Compute coordinates of the base pixel in the projection plane
    final int d0hBy4Quotient = d0h >> 2;                     assert 0 <= d0hBy4Quotient && d0hBy4Quotient <= 2;
    final int d0hMod4 = d0h & 3;                             assert 0 <= d0hMod4 && d0hMod4 <= 3;
    final int hD0h = 1 - d0hBy4Quotient;                     assert -1 <= hD0h && hD0h <= 1;
    // +1 if the base cell is not equatorial
    int lD0h = (d0hMod4 << 1) | (hD0h & 1);                  assert (hD0h == 0 && lD0h == 0 || lD0h == 2 || lD0h == 4 ||  lD0h == 6)
    || (hD0h != 0 && lD0h == 1 || lD0h == 3 || lD0h == 5 ||  lD0h == 7);
    // Compute let's go
    double lon = lInD0h * this.oneOverNside;                 assert -1 < lon && lon < 1;
    double lat = hInD0h * this.oneOverNside;                 assert -1 < lat && lat < 1;
    lon += vertexDirection.timeXOffset(this.oneOverNside);   assert -1 <= lon && lon <= 1; // The only diff with center
    lat += vertexDirection.timeYOffset(this.oneOverNside);   assert -1 <= lat && lat <= 1; // The only diff with center
    if ((d0h < 4 && hInD0h > 0) || (d0h > 7 && hInD0h < 0)) { // Polar cap
      long signLat = toBits(lat);                            assert lat == fromBits(signLat);
      lat = fromBits(signLat & BUT_SIGN_BIT_MASK_L);         assert  0 < lat && lat < 1; 
      signLat &= SIGN_BIT_MASK_L;                            assert signLat == 0 || signLat == SIGN_BIT_MASK_L;
      lat = 1 - lat;                                         assert   0 < lat && lat <  1;
      lon /= lat;                                            assert  -1 < lon && lon <  1;
      lat = fromBits(signLat | toBits(2 * acos(lat * ONE_OVER_SQRT6) - HALF_PI));
    } else { // Equatorial zone
      lat = asin((lat + hD0h) * TRANSITION_Z);
      lD0h |= ((lInD0h & SIGN_BIT_MASK_I) >>> (24 + d0h)) & 8;  assert 0 <= lD0h && lD0h <= 8;
    }
    resultLonLat[LON_INDEX] = (lon + lD0h) * PI_OVER_FOUR;
    resultLonLat[LAT_INDEX] = lat;
  }

  @Override
  public EnumMap<Cardinal, double[]> vertices(long hash, EnumSet<Cardinal> cardinalPoints) {
    final EnumMap<Cardinal, double[]> verticesMap = new EnumMap<Cardinal, double[]>(Cardinal.class);
    for (final Cardinal c : cardinalPoints) {
      final double[] resultLonLat = new double[2];
      verticesMap.put(c, resultLonLat);
    }
    this.vertices(hash, verticesMap);
    return verticesMap;
  }

  @Override
  public void vertices(long hash, EnumMap<Cardinal, double[]> cardinalPoints) {
    checkHashRange(hash);
    // Pull apart the hash elements
    final int d0h = (int) (hash >> this.twiceDepth);         assert 0 <= d0h && d0h < 12;
    hash = this.fc.hash2ij(hash & this.xyMask);
    final int iInD0h = this.fc.ij2i(hash);                   assert 0 <= iInD0h && iInD0h < this.nside;
    final int jInD0h = this.fc.ij2j(hash);                   assert 0 <= jInD0h && jInD0h < this.nside;
    // Compute coordinates from the center of the base pixel with x-axis = W-->E, y-axis = S-->N
    final int lInD0h = iInD0h - jInD0h;                      assert -nside < lInD0h && lInD0h < nside;
    final int hInD0h = iInD0h + jInD0h - this.modNsideMask;  assert -nside < hInD0h && hInD0h < nside;
    // Compute coordinates of the base pixel in the projection plane
    final int d0hBy4Quotient = d0h >> 2;                     assert 0 <= d0hBy4Quotient && d0hBy4Quotient <= 2;
    final int d0hMod4 = d0h & 3;                             assert 0 <= d0hMod4 && d0hMod4 <= 3;
    final int hD0h = 1 - d0hBy4Quotient;                     assert -1 <= hD0h && hD0h <= 1;
    // +1 if the base cell is not equatorial
    int lD0h = (d0hMod4 << 1) | (hD0h & 1);                  assert (hD0h == 0 && lD0h == 0 || lD0h == 2 || lD0h == 4 ||  lD0h == 6)
    || (hD0h != 0 && lD0h == 1 || lD0h == 3 || lD0h == 5 ||  lD0h == 7);
    // Compute let's go
    double lon = lInD0h * this.oneOverNside;                 assert -1 < lon && lon < 1;
    double lat = hInD0h * this.oneOverNside;                 assert -1 < lat && lat < 1;
    final double[] east = cardinalPoints.get(Cardinal.E);
    final double[] west = cardinalPoints.get(Cardinal.W);
    final double[] north = cardinalPoints.get(Cardinal.N);
    final double[] south = cardinalPoints.get(Cardinal.S);
    if ((d0h < 4 && hInD0h > 0) || (d0h > 7 && hInD0h < 0)) { // Polar cap
      // West-East if required
      if (east != null || west != null) {
        double latEW = lat;
        long signLat = toBits(latEW);
        latEW = fromBits(signLat & BUT_SIGN_BIT_MASK_L);
        signLat &= SIGN_BIT_MASK_L;
        final double c = 1 - latEW;
        latEW = 2 * acos(c * ONE_OVER_SQRT6) - HALF_PI;
        latEW = fromBits(signLat | toBits(latEW));
        if (east != null) {
          east[0] = ((lon + this.oneOverNside) / c + lD0h) * PI_OVER_FOUR;
          east[1] = latEW;
        }
        if (west != null) {
          west[0] = ((lon - this.oneOverNside) / c + lD0h) * PI_OVER_FOUR;
          west[1] = latEW;
        }
      }
      // North-South if required
      if (north != null) {
        double latN = lat + this.oneOverNside;
        double lonN = lon;
        long signLat = toBits(latN);
        latN = fromBits(signLat & BUT_SIGN_BIT_MASK_L);
        signLat &= SIGN_BIT_MASK_L;
        latN = 1 - latN;
        lonN /= latN;
        north[0] = (lonN + lD0h) * PI_OVER_FOUR;
        north[1] = fromBits(signLat | toBits(2 * acos(latN * ONE_OVER_SQRT6) - HALF_PI));
      } 
      if (south != null) {
        double latS = lat - this.oneOverNside;
        double lonS = lon;
        long signLat = toBits(latS);
        latS = fromBits(signLat & BUT_SIGN_BIT_MASK_L);
        signLat &= SIGN_BIT_MASK_L;
        latS = 1 - latS;
        lonS /= latS;
        south[0] = (lonS + lD0h) * PI_OVER_FOUR;
        south[1] = fromBits(signLat | toBits(2 * acos(latS * ONE_OVER_SQRT6) - HALF_PI));
      }
    } else { // Equatorial zone
      lD0h |= ((lInD0h & SIGN_BIT_MASK_I) >>> (24 + d0h)) & 8;  assert 0 <= lD0h && lD0h <= 8;
      // West-East if required
      if (east != null || west != null) {
        final double latEW = asin((lat + hD0h) * TRANSITION_Z);
        if (east != null) {
          east[0] = (lon + lD0h + this.oneOverNside) * PI_OVER_FOUR;
          east[1] = latEW;
        }
        if (west != null) {
          west[0] = (lon + lD0h - this.oneOverNside) * PI_OVER_FOUR;
          west[1] = latEW;
        }
      }
      // North-South if required
      if (north != null || south != null) {
        lon += lD0h;
        lon *= PI_OVER_FOUR;
        if (north != null) {
          north[0] = lon;
          north[1] = asin((lat + hD0h + this.oneOverNside) * TRANSITION_Z);
        }
        if (south != null) {
          south[0] = lon;
          south[1] = asin((lat + hD0h - this.oneOverNside) * TRANSITION_Z);
        }
      }
    }
  }

  @Override
  public double[][] pathAlongCellSide(long hash, Cardinal fromVertex, Cardinal toVertex,
      boolean isToVertexIncluded, int nSegments) {
    final int resultSize = isToVertexIncluded ? nSegments + 1 : nSegments;
    final double[][] pathPoints = new double[resultSize][];
    pathAlongCellSide(hash, fromVertex, toVertex, isToVertexIncluded, nSegments, pathPoints);
    return pathPoints;
  }

  @Override
  public void pathAlongCellSide(long hash, Cardinal fromVertex, Cardinal toVertex,
      boolean isToVertexIncluded, int nSegments, double[][] pathPoints) {
    pathAlongCellSide(hash, fromVertex, toVertex, isToVertexIncluded, nSegments, pathPoints, 0);
  }

  @Override
  public double[][] pathAlongCellEdge(long hash, Cardinal startingVertex,
      boolean clockwiseDirection, int nSegmentsBySide) {
    final double[][] pathPoints = new double[nSegmentsBySide << 2][];
    pathAlongCellEdge(hash, startingVertex, clockwiseDirection, nSegmentsBySide, pathPoints);
    return pathPoints;
  }

  @Override
  public void pathAlongCellEdge(long hash, Cardinal startingVertex, boolean clockwiseDirection,
      int nSegmentsBySide, double[][] pathPoints) {
    checkHashRange(hash);
    // Pull apart the hash elements
    final int d0h = (int) (hash >> this.twiceDepth);
    hash = this.fc.hash2ij(hash & this.xyMask);
    final int iInD0h = this.fc.ij2i(hash);
    final int jInD0h = this.fc.ij2j(hash);
    // Compute coordinates from the center of the base pixel with x-axis = W-->E, y-axis = S-->N
    final int lInD0h = iInD0h - jInD0h;
    final int hInD0h = iInD0h + jInD0h - this.modNsideMask;
    // Compute coordinates of the base pixel in the projection plane
    final int d0hBy4Quotient = d0h >> 2;
      final int d0hMod4 = d0h & 3;
      final int hD0h = 1 - d0hBy4Quotient;
      // +1 if the base cell is not equatorial
      int lD0h = (d0hMod4 << 1) | (hD0h & 1);
      // Compute let's go
      double lonC = lInD0h * this.oneOverNside;
      double latC = hInD0h * this.oneOverNside;
      // Unrolled loop over successive sides
      // - compute vertex sequence
      Cardinal vertex1 = startingVertex, vertex2, vertex3, vertex4;
      if (clockwiseDirection) {
        vertex2 = startingVertex.nextClockwise();
        vertex3  = vertex2.nextClockwise();
        vertex4 = vertex3.nextClockwise();
      } else {
        vertex2 = startingVertex.nextCounterClockwise();
        vertex3  = vertex2.nextCounterClockwise();
        vertex4 = vertex3.nextCounterClockwise();
      }
      // - make the five sides
      if ((d0h < 4 && hInD0h > 0) || (d0h > 7 && hInD0h < 0)) { // Polar cap
        pathAlongCellSidePolarCap(lonC, latC, lD0h, vertex1, vertex2, false, nSegmentsBySide, pathPoints, 0);
        pathAlongCellSidePolarCap(lonC, latC, lD0h, vertex2, vertex3, false, nSegmentsBySide, pathPoints, nSegmentsBySide);
        pathAlongCellSidePolarCap(lonC, latC, lD0h, vertex3, vertex4, false, nSegmentsBySide, pathPoints, nSegmentsBySide << 1);
        pathAlongCellSidePolarCap(lonC, latC, lD0h, vertex4, vertex1, false, nSegmentsBySide, pathPoints, (nSegmentsBySide << 2) - nSegmentsBySide);
      } else {
        pathAlongCellSideEquatRegion(lonC, latC, lD0h, d0h, hD0h, lInD0h, vertex1, vertex2, false, nSegmentsBySide, pathPoints, 0);
        pathAlongCellSideEquatRegion(lonC, latC, lD0h, d0h, hD0h, lInD0h, vertex2, vertex3, false, nSegmentsBySide, pathPoints, nSegmentsBySide);
        pathAlongCellSideEquatRegion(lonC, latC, lD0h, d0h, hD0h, lInD0h, vertex3, vertex4, false, nSegmentsBySide, pathPoints, nSegmentsBySide << 1);
        pathAlongCellSideEquatRegion(lonC, latC, lD0h, d0h, hD0h, lInD0h, vertex4, vertex1, false, nSegmentsBySide, pathPoints, (nSegmentsBySide << 2) - nSegmentsBySide);
      }
  }

  private void pathAlongCellSide(long hash,
      final Cardinal fromVertex, final Cardinal toVertex, 
      boolean isToVertexIncluded, int nSegments, final double[][] pathPoints, 
      final int fromPathPointsIndex) {
    checkHashRange(hash);
    // Pull apart the hash elements
    final int d0h = (int) (hash >> this.twiceDepth);
    hash = this.fc.hash2ij(hash & this.xyMask);
    final int iInD0h = this.fc.ij2i(hash);
    final int jInD0h = this.fc.ij2j(hash);
    // Compute coordinates from the center of the base pixel with x-axis = W-->E, y-axis = S-->N
    final int lInD0h = iInD0h - jInD0h;
    final int hInD0h = iInD0h + jInD0h - this.modNsideMask;
    // Compute coordinates of the base pixel in the projection plane
    final int d0hBy4Quotient = d0h >> 2;
        final int d0hMod4 = d0h & 3;
        final int hD0h = 1 - d0hBy4Quotient;
        // +1 if the base cell is not equatorial
        int lD0h = (d0hMod4 << 1) | (hD0h & 1);
        // Compute let's go
        double lonC = lInD0h * this.oneOverNside;
        double latC = hInD0h * this.oneOverNside;
        // Let's go
        if ((d0h < 4 && hInD0h > 0) || (d0h > 7 && hInD0h < 0)) { // Polar cap
          pathAlongCellSidePolarCap(lonC, latC, lD0h, fromVertex, toVertex,
              isToVertexIncluded, nSegments, pathPoints, fromPathPointsIndex);
        } else { // Equatorial zone
          pathAlongCellSideEquatRegion(lonC, latC, lD0h, d0h, hD0h, lInD0h, fromVertex, toVertex,
              isToVertexIncluded, nSegments, pathPoints, fromPathPointsIndex);
        }
  }

  private void pathAlongCellSidePolarCap(final double lonC, final double latC, final int lD0h,
      final Cardinal fromVertex, final Cardinal toVertex, 
      boolean isToVertexIncluded, int nSegments, final double[][] pathPoints, 
      final int fromPathPointsIndex) {
    final int resultSize = isToVertexIncluded ? nSegments + 1 : nSegments;
    // Compute starting point offsets
    final double fromOffsetX = fromVertex.timeXOffset(this.oneOverNside);
    final double fromOffsetY = fromVertex.timeYOffset(this.oneOverNside);
    // Compute stepX and stepY
    final double stepX = (toVertex.timeXOffset(this.oneOverNside) - fromOffsetX) / nSegments;
    final double stepY = (toVertex.timeYOffset(this.oneOverNside) - fromOffsetY) / nSegments;
    for (int i = 0; i < resultSize; i++) {
      final double[] pathPoint = new double[2];
      double lon = lonC + fromOffsetX + i * stepX;
      double lat = latC + fromOffsetY + i * stepY;
      long signLat = toBits(lat);
      lat = fromBits(signLat & BUT_SIGN_BIT_MASK_L);
      signLat &= SIGN_BIT_MASK_L;
      lat = 1 - lat;
      lon /= lat;
      lat = fromBits(signLat | toBits(2 * acos(lat * ONE_OVER_SQRT6) - HALF_PI));
      pathPoint[LON_INDEX] = (lon + lD0h) * PI_OVER_FOUR;
      pathPoint[LAT_INDEX] = lat;
      pathPoints[i + fromPathPointsIndex] = pathPoint;
    }
  }

  private void pathAlongCellSideEquatRegion(final double lonC, final double latC, int lD0h,
      final int d0h, final int hD0h, final int lInD0h,
      final Cardinal fromVertex, final Cardinal toVertex, 
      boolean isToVertexIncluded, int nSegments, final double[][] pathPoints, 
      final int fromPathPointsIndex) {
    final int resultSize = isToVertexIncluded ? nSegments + 1 : nSegments;
    // Compute starting point offsets
    final double fromOffsetX = fromVertex.timeXOffset(this.oneOverNside);
    final double fromOffsetY = fromVertex.timeYOffset(this.oneOverNside);
    // Compute stepX and stepY
    final double stepX = (toVertex.timeXOffset(this.oneOverNside) - fromOffsetX) / nSegments;
    final double stepY = (toVertex.timeYOffset(this.oneOverNside) - fromOffsetY) / nSegments;
    lD0h |= ((lInD0h & SIGN_BIT_MASK_I) >>> (24 + d0h)) & 8;
    for (int i = 0; i < resultSize; i++) {
      final double[] pathPoint = new double[2];
      double lon = lonC + fromOffsetX + i * stepX;
      double lat = latC + fromOffsetY + i * stepY;
      lat = asin((lat + hD0h) * TRANSITION_Z);
      pathPoint[LON_INDEX] = (lon + lD0h) * PI_OVER_FOUR;
      pathPoint[LAT_INDEX] = lat;
      pathPoints[i + fromPathPointsIndex] = pathPoint;
    }
  }

  ///////////////////////////////
  // Neighbours implementation //
  ///////////////////////////////

  private static byte[][] createNorhtPolarCapNeighbours(final int d0h) {
    assert 0 <= d0h && d0h < 4 : "d0h: " + d0h; // d0h = 0, 1, 2 or 3
    // Compute the offset with respect to the hash=
    final byte offset = (byte) d0h;
    assert 0 <= offset && offset <= 3 : "Offset: " + offset;
    // Empty value: all bits set to 1
    final byte x = -1;
    return createNeighbours(
        offset,                          // k = 0, C
        (byte) (  4 + offset)          , // k = 1, SW
        (byte) ( (1 + offset) & 3)     , // k = 2, NE, (y, nside-1)
        (byte) (((5 + offset) & 7) | 4), // k = 3, SE
        (byte) (  8 + offset)          , // k = 4, S
        x,                               // k = 5, E
        (byte) ( (3 + offset) & 3)     , // k = 6, NW, (nside-1, x)
        x,                               // k = 7, W
        (byte) ( (2 + offset) & 3)       // k = 8, N, (nside-1, nside-1) instead of (0, 0)
        );
  }

  private static final byte[][] createEquatorialRegionNeighbours(final int d0h) {
    assert 4 <= d0h && d0h <= 7 : "d0h: " + d0h; // d0h = 4, 5, 6 or 7
    // Compute the offset with respect to the hash=4
    final byte offset = (byte) (d0h - 4);
    assert 0 <= offset && offset <= 3 : "Offset: " + offset;
    // Empty value: all bits set to 1
    final byte x = -1;
    return createNeighbours(
        (byte) d0h,                      // k = 0, C
        // sw = 7 + (offset == 0 ? 4 : offset);
        (byte) ((11 + offset) & 11),     // k = 1, SW
        offset,                          // k = 2, NE
        (byte) (8 + offset),             // k = 3, SE
        x,                               // k = 4, S
        (byte) (((5 + offset) & 7) | 4), // k = 5, E
        (byte) ( (3 + offset) & 3)     , // k = 6, NW
        (byte) (((7 + offset) & 7) | 4), // k = 7, W
        x                                // k = 8, N
        );
  }

  private static byte[][] createSouthPolarCapNeighbours(final int d0h) {
    assert 8 <= d0h && d0h < 12 : "d0h: " + d0h; // d0h = 8, 9, 10 or 11
    // Compute the offset with respect to the hash
    final byte offset = (byte) (d0h - 8);
    assert 0 <= offset && offset <= 3 : "Offset: " + offset;
    // Empty value: all bits set to 1
    final byte x = -1;
    return createNeighbours(
        (byte) d0h,                       // k = 0, C
        (byte) ((11 + offset) & 11),      // k = 1, SW, (nside-1, x)
        (byte) (((5 + offset) &  7) | 4), // k = 2, NE
        (byte) (( 9 + offset) & 11),      // k = 3, SE, (y, nside-1)
        (byte) ((10 + offset) & 11),      // k = 4, S, (nside-1, nside-1)
        x,                                // k = 5, E
        (byte) (4 + offset),              // k = 6, NW
        x,                                // k = 7, W
        offset                            // k = 8, N
        );
  }

  private static final byte[][] createNeighbours(
      final byte  c, final byte sw, final byte ne,
      final byte se, final byte  s, final byte  e,
      final byte nw, final byte  w, final byte  n) {
    // The order in each array is:
    //          0        1        2         3       4       5        6        7
    //          s       se        e        sw      ne       w       nw        n
    // NPC:                                    (y, nside-1)     (nside-1,x)  (nside-1, nside-1)
    // EQR: (x-1,y-1) (x,y-1) (x+1,y-1) (x-1,y) (x+1,y) (x-1,y+1) (x,y+1) (x+1,y+1)
    // SPC: 
    // With NPC: North polar Cap; EQR: Equatorial Region; SPC: South polar Cap
    return new byte[][]{
      {  s, se,  e, sw, ne,  w, nw,  n}, // k=0, C (used only at depth 0)
      { sw,  c,  c, sw,  c, sw,  c,  c}, // k=1, SW
      {  c,  c, ne,  c, ne,  c,  c, ne}, // k=2, NE
      { se, se, se,  c,  c,  c,  c,  c}, // k=3, SE
      {  s, se, se, sw,  c, sw,  c,  c}, // k=4,  S
      { se, se,  e,  c, ne,  c,  c, ne}, // k=5,  E
      {  c,  c,  c,  c,  c, nw, nw, nw}, // k=6, NW
      { sw,  c,  c, sw,  c,  w, nw, nw}, // k=7,  W
      {  c,  c, ne,  c, ne, nw, nw,  n}  // k=8,  N
    };
  }

  private static final byte[] getNeighbours(final int d0h, final int k) {
    // WE DECIDED TO MAKE A STATIC VERSION, USING BYTE INSTEAD OF LONG, WITHOUT PRE << twiceDepth
    // IN ORDER TO MINIMIZE THE CPU CACHE OCCUPATION
    if (NEIGHBOURS[d0h] == null) {
      synchronized (NEIGHBOURS) {
        if (NEIGHBOURS[d0h] == null) {
          NEIGHBOURS[d0h] = d0h < 4 ? createNorhtPolarCapNeighbours(d0h)
              : d0h < 8 ? createEquatorialRegionNeighbours(d0h)
                  : createSouthPolarCapNeighbours(d0h);
        }
      }
    }
    return NEIGHBOURS[d0h][k];
  }

  
  private static void copyD0hNeig(final int d0h, final FlatHashList neighbours) {
    final byte[] a = D0H_NEIGHBOURS[d0h];
    // 4,8 or 5, 7 are null
    neighbours.hList[0] = a[0];
    neighbours.hList[1] = a[1];
    neighbours.hList[2] = a[2];
    neighbours.hList[3] = a[3];
    neighbours.hList[4] = a[4];
    neighbours.hList[5] = a[5];
    neighbours.size = 6;
  }
  
  @Override
  public void neighbours(long hash, final FlatHashList neighbours) {    
    // Compute the value of the depth 0 pixel.
    final int d0h = (int) (hash >>> this.twiceDepth);     assert 0 <= d0h && d0h < 12;
    if (this.depth == 0) {
      copyD0hNeig(d0h, neighbours);
      return;
    }
    
    // Separate the bits coding the value of the base pixel, of x, and of y
    final long d0hFull = hash & this.d0Mask;
    final long xFull   = hash & this.xMask;
    final long yFull   = hash & this.yMask;
    
    // Unshuffle to obtain xy
    hash = this.fc.hash2ij(hash & this.xyMask);
    int x = this.fc.ij2i(hash);   assert 0 <= x && x < this.nside : "x: " + x;
    int y = this.fc.ij2j(hash);   assert 0 <= y && y < this.nside : "y: " + y;

    // Define i such as i = (x == 0) ? 1 : (x == nside - 1) ? 2 : 0; 
    int i = (x - 1) >>> 31;            // i = (x==0) ? 1 : 0
    //System.out.println("d:" + depth + "; x: " + x + "; i: " + i);
    i |= ((x + 1) >> this.depth) << 1; // i |= (x==nside -1) ? 2 : 0
    //System.out.println("d:" + depth + "; x: " + x + "; i: " + i);

    assert i == 0 || i == 1 || i == 2 : "i: " + i;
    // Define j such as j = (y == 0) ? 1 : (y == nside - 1) ? 2 : 0;
    int j = (y - 1) >>> 31;
    j |= ((y + 1) >> this.depth) << 1;
    assert j == 0 || j == 1 || j == 2 : "j: " + j;

    hash = this.fc.ij2hash((x - 1) & this.modNsideMask, (y - 1) & this.modNsideMask);
    final long xm1Full = hash & this.xMask;
    final long ym1Full = hash & this.yMask;
    hash = this.fc.ij2hash((x + 1) & this.modNsideMask, (y + 1) & this.modNsideMask);
    final long xp1Full = hash & this.xMask;
    final long yp1Full = hash & this.yMask;

    if (i == 0 && j == 0) { // Not on edge, not on corner, depth > 0
      neighbours.hList[0] = d0hFull | ym1Full | xm1Full;
      neighbours.hList[1] = d0hFull | ym1Full | xFull;
      neighbours.hList[2] = d0hFull | ym1Full | xp1Full;
      neighbours.hList[3] = d0hFull | yFull   | xm1Full;
      neighbours.hList[4] = d0hFull | yFull   | xp1Full;
      neighbours.hList[5] = d0hFull | yp1Full | xm1Full;
      neighbours.hList[6] = d0hFull | yp1Full | xFull;
      neighbours.hList[7] = d0hFull | yp1Full | xp1Full;
      neighbours.size = 8;
    } else if (d0h < 4) { // d0h = 0, 1, 2 or 3, i.e. North polar cap
      final int k = (j << 1) + j + i; // <=> 3 * j + i
      assert 0 < k && k < 9 : "k: " + k; 

      final byte[] d0hneigs = getNeighbours(d0h, k);      
      int l = 0;
      switch(k) {
      case 1: case  3: case 4: // SW, SE, S
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | yFull   | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | yFull   | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | yp1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | yp1Full | xp1Full;
        break;
      case 2: // NE
        assert x == this.nside - 1;
        assert xFull == this.xMask && xp1Full == 0;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | this.yMask | (ym1Full >>> 1);
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | yFull   | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | this.yMask | (yFull >>> 1);
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | yp1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | this.yMask | (yp1Full >> 1);
        break;
      case 5: // E
        assert x == this.nside - 1 && y == 0;
        assert xFull == this.xMask && xp1Full == 0;
        assert ym1Full == this.yMask && yFull == 0 && yp1Full == 2L;
        assert (yp1Full >> 1) == 1L;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full | xFull;
        // no neighbour for k=2, i.e. for (x+1, y-1)
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | this.yMask;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | yp1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | this.yMask | 1L;
        break;
      case 6: // NW
        assert y == this.nside -1;
        assert yFull == this.yMask && yp1Full == 0;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | yFull   | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | yFull   | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | (xm1Full << 1) | this.xMask;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | (  xFull << 1) | this.xMask;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | (xp1Full << 1) | this.xMask;
        break;
      case 7: // W
        assert x == 0 && y == this.nside - 1;
        assert xFull == 0 && xm1Full == this.xMask;
        assert (xp1Full << 1) == 2L;
        assert yFull == this.yMask && yp1Full == 0;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | yFull   | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | yFull   | xp1Full;
        // no neighbour for k=5, i.e. (x-1, y+1)
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | 2L | xm1Full;
        break;
      case 8: // N: Pixel located the north pole
        assert x == this.nside - 1 && y == this.nside - 1;
        assert xFull == this.xMask && xp1Full == 0;
        assert yFull == this.yMask && yp1Full == 0;
        assert (yFull >>> 1) == xFull && (xFull << 1) == yFull;
        assert (xm1Full << 1) == ym1Full && (ym1Full >>> 1) == xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | yFull | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | yFull | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | this.xyMask; //& -1L;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | ym1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | this.xyMask; //& -1L;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | this.xyMask; //& -1L;
        break;
      default:
        assert false;
      }
      neighbours.size = l;
    } else if (d0h > 7) { // d0h = 8, 9, 10 or 11, i.e. South polar cap
      final int k = j * 3 + i;
      assert 0 < k && k < 9 : "k: " + k;
      
      final byte[] d0hneigs = getNeighbours(d0h, k);      
      int l = 0;
      switch(k) {
      case 2: case 6: case 8: // NE, NW, N
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | yFull   | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | yFull   | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | yp1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | yp1Full | xp1Full;
        break;
      case 1: // SW
        assert x == 0 && xFull == 0L;
        assert xm1Full == this.xMask && xp1Full == 1L;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | (ym1Full >> 1);
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | (yFull >> 1);
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | yFull   | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | (yp1Full >> 1);
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | yp1Full | xp1Full;
        break;
      case 3: // SE
        assert y == 0 && yFull == 0L;
        assert ym1Full == this.yMask && yp1Full == 2L;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | (xm1Full << 1);
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | (xFull << 1);
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | (xp1Full << 1);
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | yp1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | yp1Full | xp1Full;
        break;
      case 4: // S
        assert x == 0 && y == 0 && xFull == 0L && yFull == 0L;
        assert xm1Full == this.xMask && xp1Full == 1L;
        assert ym1Full == this.yMask && yp1Full == 2L;
        assert (yp1Full >> 1) == xp1Full;
        assert (xp1Full | yp1Full) == 3L;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth; //& -1L;
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | yp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | 3L;
        break;
      case 5: // E
        assert x == this.nside - 1 && y == 0;
        assert xFull == this.xMask && xp1Full == 0;
        assert ym1Full == this.yMask && yFull == 0 && yp1Full == 2L
            :( ym1Full ==  this.yMask) + " " +  (yFull == 0) + " " +  (yp1Full == 1L);
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | (xm1Full << 1);
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | this.yMask;
        // no neighbour for k=2, i.e. for (x+1, y-1)
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth;
        neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | yp1Full | xm1Full;
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full | xFull;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | yp1Full ;
        break;
      case 7: // W
        assert x == 0 && y == this.nside - 1;
        assert xFull == 0 && xm1Full == this.xMask;
        assert yFull == this.yMask && yp1Full == 0;
        assert (xp1Full << 1) == 2L;
        neighbours.hList[l++] = ((long) d0hneigs[0]) << this.twiceDepth | (ym1Full >> 1);
        neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full;
        neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
        neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;//this.xyMask; //| -1L;
        neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | yFull   | xp1Full;
        // no neighbour for k=5, i.e. (x-1, y+1)
        neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth;
        neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | xp1Full;
        break;
      default:
        assert false;
      }
      neighbours.size = l;
    } else { // d0h = 4, 5, 6 or 7, i.e. Equatorial region
      final int k = j * 3 + i;
      assert 0 < k && k < 9 : "k: " + k; 
      final byte[] d0hneigs = getNeighbours(d0h, k);
      neighbours.hList[0]   = ((long) d0hneigs[0]) << this.twiceDepth | ym1Full | xm1Full;
      // if k=4, no value associated to (x-1, y-1), so we overwrite it
      int l = (k == 4) ? 0 : 1; // <=> ((((k + 3) & 7) + 1) >> 3);
      neighbours.hList[l++] = ((long) d0hneigs[1]) << this.twiceDepth | ym1Full | xFull;
      neighbours.hList[l++] = ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
      neighbours.hList[l++] = ((long) d0hneigs[3]) << this.twiceDepth | yFull   | xm1Full;
      neighbours.hList[l++] = ((long) d0hneigs[4]) << this.twiceDepth | yFull   | xp1Full;
      neighbours.hList[l++] = ((long) d0hneigs[5]) << this.twiceDepth | yp1Full | xm1Full;
      neighbours.hList[l++] = ((long) d0hneigs[6]) << this.twiceDepth | yp1Full | xFull;
      neighbours.hList[l++] = ((long) d0hneigs[7]) << this.twiceDepth | yp1Full | xp1Full;
      // if k=8, no value associated to (x+1, y+1), so we ignore it
      l -= (k >> 3); // l -= (k == 8) ? 1 : 0;
      neighbours.size = l;
    }
  }

  @Override
  public long neighbour(long hash, MainWind direction) {
    // Compute the value of the depth 0 pixel.
    final int d0h = (int) (hash >>> this.twiceDepth);     assert 0 <= d0h && d0h < 12;
    if (this.depth == 0) {
      final byte[] d0hNeig = getNeighbours(d0h, 0);
      switch(direction) {
      case  S: return d0hNeig[0];
      case SE: return d0hNeig[1];
      case  E: return d0hNeig[2];
      case SW: return d0hNeig[3];
      case NE: return d0hNeig[4];
      case  W: return d0hNeig[5];
      case NW: return d0hNeig[6];
      case  N: return d0hNeig[7];
      default: throw new IllegalArgumentException("Wrong direction: " + direction);
      }
    }
 
    // Separate the bits coding the value of the base pixel, of x, and of y
    final long d0hFull = hash & this.d0Mask;
    final long xFull   = hash & this.xMask;
    final long yFull   = hash & this.yMask;
    
    // Unshuffle to obtain xy
    hash = this.fc.hash2ij(hash & this.xyMask);
    int x = this.fc.ij2i(hash);    assert 0 <= x && x < this.nside;
    int y = this.fc.ij2j(hash);   assert 0 <= y && y < this.nside;

    // Define i such as i = (x == 0) ? 1 : (x == nside - 1) ? 2 : 0;  
    int i = (x - 1) >>> 31;            // i = (x==0) ? 1 : 0
    i |= ((x + 1) >> this.depth) << 1; // i |= (x==nside -1) ? 2 : 0
    assert i == 0 || i == 1 || i == 2 : "i: " + i;
    // Define j such as j = (y == 0) ? 1 : (y == nside - 1) ? 2 : 0;
    int j = (y - 1) >>> 31;
    j |= ((y + 1) >> this.depth) << 1;
    assert j == 0 || j == 1 || j == 2 : "j: " + j;

    hash = this.fc.ij2hash((x - 1) & this.modNsideMask, (y - 1) & this.modNsideMask);
    final long xm1Full = hash & this.xMask;
    final long ym1Full = hash & this.yMask;
    hash = this.fc.ij2hash((x + 1) & this.modNsideMask, (y + 1) & this.modNsideMask);
    final long xp1Full = hash & this.xMask;
    final long yp1Full = hash & this.yMask;

    if (i == 0 && j == 0) { // Not on edge, not on corner, depth > 0
      switch(direction) {
      case  S: return bits2hash(d0hFull, ym1Full, xm1Full);
      case SE: return bits2hash(d0hFull, ym1Full, xFull);
      case  E: return bits2hash(d0hFull, ym1Full, xp1Full);
      case SW: return bits2hash(d0hFull, yFull  , xm1Full);
      case NE: return bits2hash(d0hFull, yFull  , xp1Full);
      case  W: return bits2hash(d0hFull, yp1Full, xm1Full);
      case NW: return bits2hash(d0hFull, yp1Full, xFull);
      case  N: return bits2hash(d0hFull, yp1Full, xp1Full);
      default: throw new IllegalArgumentException("Wrong direction: " + direction);
      }
    } else if (d0h < 4) { // d0h = 0, 1, 2 or 3, i.e. North polar cap
      final int k = (j << 1) + j + i; // <=> 3 * j + i
      assert 0 < k && k < 9 : "k: " + k; 
      final byte[] d0hneigs = getNeighbours(d0h, k);
      switch(k) {
      case 1: case  3: case 4: // SW, SE, S
        switch(direction) {
        case  S: return bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full, xm1Full);
        case SE: return bits2hash(((long) d0hneigs[1]) << this.twiceDepth, ym1Full, xFull);
        case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, ym1Full, xp1Full);
        case SW: return bits2hash(((long) d0hneigs[3]) << this.twiceDepth, yFull  , xm1Full);
        case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, yFull  , xp1Full);
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, yp1Full, xm1Full);
        case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth, yp1Full, xFull);
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, yp1Full, xp1Full);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 2: // NE
        assert x == this.nside - 1;
        assert xFull == this.xMask && xp1Full == 0;
        switch(direction) {
        case  S: return bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full   , xm1Full);
        case SE: return bits2hash(((long) d0hneigs[1]) << this.twiceDepth, ym1Full   , xFull);
        case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, this.yMask, ym1Full >>> 1);
        case SW: return bits2hash(((long) d0hneigs[3]) << this.twiceDepth, yFull     , xm1Full);
        case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, this.yMask, yFull >>> 1);
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, yp1Full   , xm1Full);
        case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth, yp1Full   , xFull);
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, this.yMask, yp1Full >> 1);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 5: // E
        assert x == this.nside - 1 && y == 0;
        assert xFull == this.xMask && xp1Full == 0;
        assert ym1Full == this.yMask && yFull == 0 && yp1Full == 2L;
        assert (yp1Full >> 1) == 1L;
        switch(direction) {
        case  S: return bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full   , xm1Full);
        case SE: return bits2hash(((long) d0hneigs[1]) << this.twiceDepth, ym1Full   , xFull);
        case  E: return -1L;
        case SW: return ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;
        case NE: return ((long) d0hneigs[4]) << this.twiceDepth | this.yMask;
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, yp1Full   , xm1Full);
        case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth, yp1Full   , xFull);
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, this.yMask, 1L);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 6: // NW
        assert y == this.nside -1;
        assert yFull == this.yMask && yp1Full == 0;
        switch(direction) {
        case  S: return bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full     , xm1Full);
        case SE: return bits2hash(((long) d0hneigs[1]) << this.twiceDepth, ym1Full     , xFull);
        case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, ym1Full     , xp1Full);
        case SW: return bits2hash(((long) d0hneigs[3]) << this.twiceDepth, yFull       , xm1Full);
        case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, yFull       , xp1Full);
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, xm1Full << 1, this.xMask);
        case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth,   xFull << 1, this.xMask);
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, xp1Full << 1, this.xMask);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 7: // W
        assert x == 0 && y == this.nside - 1;
        assert xFull == 0 && xm1Full == this.xMask;
        assert (xp1Full << 1) == 2L;
        assert yFull == this.yMask && yp1Full == 0;
        switch(direction) {
        case  S: return bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full     , xm1Full);
        case SE: return ((long) d0hneigs[1]) << this.twiceDepth | ym1Full;
        case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, ym1Full     , xp1Full);
        case SW: return bits2hash(((long) d0hneigs[3]) << this.twiceDepth, yFull       , xm1Full);
        case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, yFull       , xp1Full);
        case  W: return -1L;
        case NW: return ((long) d0hneigs[6]) << this.twiceDepth | xm1Full;
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, 2L, xm1Full);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 8: // N: Pixel located the north pole
        assert x == this.nside - 1 && y == this.nside - 1;
        assert xFull == this.xMask && xp1Full == 0;
        assert yFull == this.yMask && yp1Full == 0;
        assert (yFull >>> 1) == xFull && (xFull << 1) == yFull;
        assert (xm1Full << 1) == ym1Full && (ym1Full >>> 1) == xm1Full;
        switch(direction) {
        case  S: return bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full, xm1Full);
        case SE: return bits2hash(((long) d0hneigs[1]) << this.twiceDepth, ym1Full, xFull);
        case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, yFull  , xm1Full);
        case SW: return bits2hash(((long) d0hneigs[3]) << this.twiceDepth, yFull  , xm1Full);
        case NE: return ((long) d0hneigs[4]) << this.twiceDepth | this.xyMask;
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, ym1Full, xFull);
        case NW: return ((long) d0hneigs[6]) << this.twiceDepth | this.xyMask;
        case  N: return ((long) d0hneigs[7]) << this.twiceDepth | this.xyMask;
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      default:
        throw new IllegalArgumentException("Wrong k: " + k);
      }
    } else if (d0h > 7) { // d0h = 8, 9, 10 or 11, i.e. South polar cap
      final int k = j * 3 + i;
      assert 0 < k && k < 9 : "k: " + k; 
      final byte[] d0hneigs = getNeighbours(d0h, k);
      switch(k) {
      case 2: case 6: case 8: // NE, NW, N
        switch(direction) {
        case  S: return bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full, xm1Full);
        case SE: return bits2hash(((long) d0hneigs[1]) << this.twiceDepth, ym1Full, xFull);
        case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, ym1Full, xp1Full);
        case SW: return bits2hash(((long) d0hneigs[3]) << this.twiceDepth, yFull  , xm1Full);
        case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, yFull  , xp1Full);
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, yp1Full, xm1Full);
        case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth, yp1Full, xFull);
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, yp1Full, xp1Full);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 1: // SW
        assert x == 0 && xFull == 0L;
        assert xm1Full == this.xMask && xp1Full == 1L;
        switch(direction) {
        case  S: return ((long) d0hneigs[0]) << this.twiceDepth | (ym1Full >> 1);
        case SE: return ((long) d0hneigs[1]) << this.twiceDepth | ym1Full;
        case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, ym1Full, xp1Full);
        case SW: return ((long) d0hneigs[3]) << this.twiceDepth | (yFull >> 1);
        case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, yFull  , xp1Full);
        case  W: return ((long) d0hneigs[5]) << this.twiceDepth | (yp1Full >> 1);
        case NW: return ((long) d0hneigs[6]) << this.twiceDepth | yp1Full;
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, yp1Full, xp1Full);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 3: // SE
        assert y == 0 && yFull == 0L;
        assert ym1Full == this.yMask && yp1Full == 2L;
        switch(direction) {
        case  S: return ((long) d0hneigs[0]) << this.twiceDepth | (xm1Full << 1);
        case SE: return ((long) d0hneigs[1]) << this.twiceDepth | (xFull << 1);
        case  E: return ((long) d0hneigs[2]) << this.twiceDepth | (xp1Full << 1);
        case SW: return ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;
        case NE: return ((long) d0hneigs[4]) << this.twiceDepth | xp1Full;
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, yp1Full, xm1Full);
        case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth, yp1Full, xFull);
        case  N: return bits2hash(((long) d0hneigs[7]) << this.twiceDepth, yp1Full, xp1Full);
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 4: // S
        assert x == 0 && y == 0 && xFull == 0L && yFull == 0L;
        assert xm1Full == this.xMask && xp1Full == 1L;
        assert ym1Full == this.yMask && yp1Full == 2L;
        assert (yp1Full >> 1) == xp1Full;
        assert (xp1Full | yp1Full) == 3L;
        switch(direction) {
        case  S: return ((long) d0hneigs[0]) << this.twiceDepth;
        case SE: return ((long) d0hneigs[1]) << this.twiceDepth;
        case  E: return ((long) d0hneigs[2]) << this.twiceDepth | yp1Full;
        case SW: return ((long) d0hneigs[3]) << this.twiceDepth;
        case NE: return ((long) d0hneigs[4]) << this.twiceDepth | xp1Full;
        case  W: return ((long) d0hneigs[5]) << this.twiceDepth | xp1Full;
        case NW: return ((long) d0hneigs[6]) << this.twiceDepth | yp1Full;
        case  N: return ((long) d0hneigs[7]) << this.twiceDepth | 3L;
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 5: // E
        assert x == this.nside - 1 && y == 0;
        assert xFull == this.xMask && xp1Full == 0;
        assert ym1Full == this.yMask && yFull == 0 && yp1Full == 2L
            :( ym1Full ==  this.yMask) + " " +  (yFull == 0) + " " +  (yp1Full == 1L);
        switch(direction) {
        case  S: return ((long) d0hneigs[0]) << this.twiceDepth | (xm1Full << 1);
        case SE: return ((long) d0hneigs[1]) << this.twiceDepth | this.yMask;
        case  E: return -1L;
        case SW: return ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;
        case NE: return ((long) d0hneigs[4]) << this.twiceDepth;
        case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, yp1Full, xm1Full);
        case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth, yp1Full, xFull);
        case  N: return ((long) d0hneigs[7]) << this.twiceDepth | yp1Full;
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      case 7: // W
        assert x == 0 && y == this.nside - 1;
        assert xFull == 0 && xm1Full == this.xMask;
        assert yFull == this.yMask && yp1Full == 0;
        assert (xp1Full << 1) == 2L;
        switch(direction) {
        case  S: return ((long) d0hneigs[0]) << this.twiceDepth | (ym1Full >> 1);
        case SE: return ((long) d0hneigs[1]) << this.twiceDepth | ym1Full;
        case  E: return ((long) d0hneigs[2]) << this.twiceDepth | ym1Full | xp1Full;
        case SW: return ((long) d0hneigs[3]) << this.twiceDepth | xm1Full;
        case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, yFull  , xp1Full);
        case  W: return -1L;
        case NW: return ((long) d0hneigs[6]) << this.twiceDepth;
        case  N: return ((long) d0hneigs[7]) << this.twiceDepth | xp1Full;
        default: throw new IllegalArgumentException("Wrong direction: " + direction);
        }
      default:
        throw new IllegalArgumentException("Wrong k: " + k);
      }
    } else { // d0h = 4, 5, 6 or 7, i.e. Equatorial region
      final int k = j * 3 + i;
      assert 0 < k && k < 9 : "k: " + k; 
      final byte[] d0hneigs = getNeighbours(d0h, k);
      switch(direction) {
      case  S: return (k == 4) ? -1L : bits2hash(((long) d0hneigs[0]) << this.twiceDepth, ym1Full, xm1Full);
      case SE: return bits2hash(((long) d0hneigs[1]) << this.twiceDepth, ym1Full, xFull);
      case  E: return bits2hash(((long) d0hneigs[2]) << this.twiceDepth, ym1Full, xp1Full);
      case SW: return bits2hash(((long) d0hneigs[3]) << this.twiceDepth, yFull  , xm1Full);
      case NE: return bits2hash(((long) d0hneigs[4]) << this.twiceDepth, yFull  , xp1Full);
      case  W: return bits2hash(((long) d0hneigs[5]) << this.twiceDepth, yp1Full, xm1Full);
      case NW: return bits2hash(((long) d0hneigs[6]) << this.twiceDepth, yp1Full, xFull);
      case  N: return (k == 8) ? -1L : bits2hash(((long) d0hneigs[7]) << this.twiceDepth, yp1Full, xp1Full);
      default: throw new IllegalArgumentException("Wrong direction: " + direction);
      }
    }
  }

  @Override
  public NeighbourList neighbours(long hash) {
    final NeighbourList result = new NeighbourList(this.depth);
    this.neighbours(hash, result);
    return result;
  }


  @Override
  public NeighbourList neighbours(long hash, EnumSet<MainWind> directions) {
    final NeighbourList result = new NeighbourList(this.depth);
    this.neighbours(hash, directions, result);
    return result;
  }
  @Override
  public void neighbours(long hash, NeighbourList result) {
    checkHashRange(hash);
    result.clear();
    // Separate the bits coding the value of the base pixel, of x, and of y
    final long d0hBits = hash & this.d0Mask;
    final long iInD0hBits   = hash & this.xMask;
    final long jInD0hBits   = hash & this.yMask;

    if (isInBasePixelBorderFromBits(iInD0hBits, jInD0hBits)) {
      edgePixelNeighbours(hash, result);
    } else {
      innerPixelNeighbours(d0hBits, iInD0hBits, jInD0hBits, result);
    }
  }
  @Override
  public void neighbours(long hash, EnumSet<MainWind> directions, NeighbourList result) {
    checkHashRange(hash);
    // Separate the bits coding the value of the base pixel, of x, and of y
    final long d0hBits = hash & this.d0Mask;
    final long iInD0hBits   = hash & this.xMask;
    final long jInD0hBits   = hash & this.yMask;

    result.clear();
    if (isInBasePixelBorderFromBits(iInD0hBits, jInD0hBits)) {
      edgePixelNeighbours(hash, directions, result);
    } else {
      innerPixelNeighbours(d0hBits, iInD0hBits, jInD0hBits, directions, result);
    }
  }

  @Override
  public FlatHashList internalEdges(long hash, int toEdgeDeltaDepth) {
    assert 1 < toEdgeDeltaDepth && toEdgeDeltaDepth < 30;
    final int n = ((1 << toEdgeDeltaDepth) - 1) << 2; // 4 * nside - 4
    final FlatHashList result = new FlatHashList(this.depth + toEdgeDeltaDepth, n);
    internalEdges(hash, toEdgeDeltaDepth, result);
    return result;
  }
  @Override
  public void internalEdges(long hash, int toEdgeDeltaDepth, FlatHashList result) {
    // Compute the x and y part masks for deltaDepth.
    final long xMaskDD = xMask(toEdgeDeltaDepth);
    final long yMaskDD = yMask(toEdgeDeltaDepth);
    // Prepare hashes of depth of (this.depth + deltaDepth), switching hash
    // bits of 2 deltaDepth to the left.
    hash <<= (toEdgeDeltaDepth << 1);
    assert (toEdgeDeltaDepth << 1) == 2 * toEdgeDeltaDepth;
    // Prepare filling the result.
    // am1 stands for a - 1, i.e. nSide - 1,
    //  i.e. the index of the last cell along the x or y-axis
    final int am1 = (1 << toEdgeDeltaDepth) - 1; // 2^deltaDepth - 1
    assert (1 << toEdgeDeltaDepth) == pow(2, toEdgeDeltaDepth);
    int k1 = am1, k2 = (am1 << 1), k3 = am1 << 2;
    // Put the values of the 4 corners
    result.hList[0] = hash;
    result.hList[k1++] = hash | xMaskDD;
    result.hList[k2] = hash | yMaskDD | xMaskDD;
    k2 += k1;
    result.hList[--k2] = hash | yMaskDD;
    // Set the 4 sides in a single for loop
    for (int k0 = 1; k0 < am1; k0++) {
      final long kx = this.fc.i02hash(k0);//.xy2lhash(k0); // we know i * i hold on an integer
      final long ky = kx << 1;
      // Southeast axis, i.e. x-axis, y = 0
      result.hList[k0] = hash | kx;
      // Notheast axis, i.e. y-axis, x = am1
      result.hList[k1++] = hash | ky | xMaskDD;
      // Northwest axis, i.e. x-axis, y = am1
      result.hList[--k2] = hash | yMaskDD | kx;
      // Southwest axis, i.e. y-axis, x = 0
      result.hList[--k3] = hash | ky;
    }
    result.size = am1 << 2;
    assert result.size == 4 * am1;
  }

  @Override
  public FlatHashList sortedInternalEdges(final long hash, final int deltaDepth) {
    assert 1 < deltaDepth && deltaDepth < 30;
    final int n = ((1 << deltaDepth) - 1) << 1; // 4 * nside - 4
    final FlatHashList result = new FlatHashList(this.depth + deltaDepth, n);
    sortedInternalEdges(hash, deltaDepth, result);
    return result;
  }

  @Override
  public void sortedInternalEdges(long hash, final int toEdgeDeltaDepth, final FlatHashList result) {
    // Compute the x and y part masks for deltaDepth.
    final long xMaskDD = xMask(toEdgeDeltaDepth);
    final long yMaskDD = yMask(toEdgeDeltaDepth);
    // Prepare depth of order this.depth + deltaDepth hashes.
    hash <<= (toEdgeDeltaDepth << 1);    assert (toEdgeDeltaDepth << 1) == 2 * toEdgeDeltaDepth;
    // Set grid size (nSide inside the cell of depth this.depth)
    final int nSide = (1 << toEdgeDeltaDepth);
    final int am1 = nSide - 1;
    final int nHalfSide = nSide >>> 1;
    // South sub-square (dividing in 4 sub-squares)
    int x = 1, tmp;
    int k0 = 0, lim = 2, k1 = 2, k2 = am1 + nHalfSide, k3 = (am1 << 1) + nHalfSide;
    int size = am1 << 2;
    // Set South corner (first element)
    result.hList[k0++] = hash;
    // Set east corner
    result.hList[k2 - 1] = hash | xMaskDD;
    // Set west corner
    result.hList[k3 - 1] = hash | yMaskDD;
    // Set north corner (last eslement)
    result.hList[size - k0] = hash | yMaskDD | xMaskDD;
    for (; x < nHalfSide;) { // while (k < nHalfSize)
      long xs = this.fc.ij2hash(x++, nSide - x); // x shuffled
      final long xn = xs & yMaskDD;
      xs &= xMaskDD;     
      // South square, south east part
      result.hList[k0++] = hash | xs;
      // South square, south west part
      result.hList[k1++] = hash | (xs << 1);
      // East square, nort east part
      result.hList[k2++] = hash | (xs << 1) | xMaskDD;
      // West square, north west part
      result.hList[k3++] = hash | yMaskDD | xs;
      // North square, north west
      result.hList[size - k0] = hash | yMaskDD | (xn >> 1);
      // North square, north EAST
      result.hList[size - k1] = hash | xn | xMaskDD;
      // West square, north west part
      result.hList[size - k2] = hash | xn;
      // East square, nort east part
      result.hList[size - k3] = hash | (xn >> 1);
      // Change k0, k1 and limit if x== limit.
      // The following lines of code are equivalent to:
      /* if (x == lim) {
              k0 = k1;
              k1 += lim; // +2 +4 +8
              lim <<= 1; // 4 8 16 32 ...
          } */
      // To be tested if they are faster (since no risk of branch miss-prediction):
      // probably true for small deltaDepth but not for large deltaDepth.
      tmp = x & lim;     assert (x < lim && tmp == 0) || (x == lim && tmp == x);
      k0 += (tmp >> 1);
      k1 += tmp;
      tmp -= x;          assert (x < lim && tmp <  0) || (x == lim && tmp == 0);
      tmp = 1 >> tmp;    assert (x < lim && tmp == 0) || (x == lim && tmp == 1);
      lim <<= tmp;
    }
    result.size = am1 << 2;
    assert result.size == 4 * (nSide - 1);
  }

  @Override
  public FlatHashList sortedInternalEdge(final long hash, final int deltaDepth,
      final Ordinal direction) {
    final FlatHashList result = new FlatHashList(this.depth + deltaDepth, 1 << deltaDepth);
    sortedInternalEdge(hash, deltaDepth, direction, result);
    return result;
  }

  @Override
  public void sortedInternalEdge(final long hash, final int deltaDepth,
      final Ordinal direction, final FlatHashList result) {
    direction.orderedInternalEdge(this, hash, deltaDepth, result);
  }

  @Override
  public FlatHashList sortedInternalEdgeSE(long hash, final int deltaDepth) {
    final FlatHashList result = new FlatHashList(this.depth + deltaDepth, 1 << deltaDepth);
    sortedInternalEdgeSE(hash, deltaDepth, result);
    return result;
  }

  @Override
  public void sortedInternalEdgeSE(long hash, final int deltaDepth,
      final FlatHashList result) {
    hash <<= (deltaDepth << 1);
    final int nSide = 1 << deltaDepth; // 2^deltaDepth
    for (int x = 0; x < nSide; x++) {
      result.hList[x] = hash | this.fc.i02hash(x);
    }
    result.size = nSide;
  }

  @Override
  public FlatHashList sortedInternalEdgeNE(long hash, final int deltaDepth) {
    final FlatHashList result = new FlatHashList(this.depth + deltaDepth, 1 << deltaDepth);
    sortedInternalEdgeNE(hash, deltaDepth, result);
    return result;
  }

  @Override
  public void sortedInternalEdgeNE(long hash, final int deltaDepth,
      final FlatHashList result) {
    hash <<= (deltaDepth << 1);
    final int nSide = 1 << deltaDepth; // 2^deltaDepth
    final long x = (nSide - 1) << 32;
    for (int y = 0; y < nSide; y++) {
      result.hList[y] = hash | this.fc.xy2hash(x >> 32, y /*y | x*/);
    }
    result.size = nSide;
  }

  @Override
  public FlatHashList sortedInternalEdgeNW(long hash, final int deltaDepth) {
    final FlatHashList result = new FlatHashList(this.depth + deltaDepth, 1 << deltaDepth);
    sortedInternalEdgeNW(hash, deltaDepth, result);
    return result;
  }

  @Override
  public void sortedInternalEdgeNW(long hash, final int deltaDepth,
      final FlatHashList result) {
    hash <<= (deltaDepth << 1);
    final int nSide = 1 << deltaDepth; // 2^deltaDepth
    final long y = (nSide - 1) << 32;
    for (int x = 0; x < nSide; x++) {
      result.hList[x] = hash | this.fc.xy2hash(x, y >> 32/*y | x*/);
    }
    result.size = nSide;
  }

  @Override
  public FlatHashList sortedInternalEdgeSW(long hash, final int deltaDepth) {
    final FlatHashList result = new FlatHashList(this.depth + deltaDepth, 1 << deltaDepth);
    sortedInternalEdgeSW(hash, deltaDepth, result);
    return result;
  }

  @Override
  public void sortedInternalEdgeSW(long hash, final int deltaDepth,
      final FlatHashList result) {
    hash <<= (deltaDepth << 1);
    final int nSide = 1 << deltaDepth; // 2^deltaDepth
    for (int y = 0; y < nSide; y++) {
      result.hList[y] = hash | this.fc.i02hash(y);
    }
    result.size = nSide;
  }

  @Override
  public long internalCorner(long hash, int toEdgeDeltaDepth, Cardinal direction) {
    return direction.internalCorner(this, hash, toEdgeDeltaDepth);
  }
  @Override
  public long internalCornerN(long hash, int toEdgeDeltaDepth) {
    hash <<= (toEdgeDeltaDepth << 1);
    return hash | xyMask(toEdgeDeltaDepth);
  }
  @Override
  public long internalCornerS(long hash, int toEdgeDeltaDepth) {
    hash <<= (toEdgeDeltaDepth << 1);
    return hash;
  }
  @Override
  public long internalCornerE(long hash, int toEdgeDeltaDepth) {
    hash <<= (toEdgeDeltaDepth << 1);
    return hash | yMask(toEdgeDeltaDepth);
  }
  @Override
  public long internalCornerW(long hash, int toEdgeDeltaDepth) {
    hash <<= (toEdgeDeltaDepth << 1);
    return hash | xMask(toEdgeDeltaDepth);
  }

  @Override
  public FlatHashList externalEdges(long hash, int deltaDepth) {
    final FlatHashList result =  new FlatHashList(this.depth + deltaDepth, 4 + (4 << deltaDepth));
    externalEdges(hash, deltaDepth, result);
    return result;
  }

  @Override
  public void externalEdges(long hash, int toEdgeDeltaDepth, final FlatHashList result) {
    externalEdges(hash, toEdgeDeltaDepth, result, false);
  }

  @Override
  public FlatHashList sortedExternalEdges(long hash, int deltaDepth) {
    final FlatHashList result =  new FlatHashList(this.depth + deltaDepth, 4 + (4 << deltaDepth));
    sortedExternalEdges(hash, deltaDepth, result);
    return result;
  }

  @Override
  public void sortedExternalEdges(long hash, int toEdgeDeltaDepth, final FlatHashList result) {
    externalEdges(hash, toEdgeDeltaDepth, result, true);
  }

  private void externalEdges(long hash, int toEdgeDeltaDepth, final FlatHashList result, boolean sorted) {
    checkHashRange(hash);
    result.clear();
    final NeighbourList neihbours = new NeighbourList(this.depth);
    // Separate the bits coding the value of the base pixel, of x, and of y
    final long d0hBits = hash & this.d0Mask;
    final long iInD0hBits   = hash & this.xMask;
    final long jInD0hBits   = hash & this.yMask;

    if (isInBasePixelBorderFromBits(iInD0hBits, jInD0hBits)) {
      // Not easy: opposite directions depends on base cell neighbours
      final int d0h = (int) (hash >> this.twiceDepth);
      hash = this.fc.hash2ij(hash & this.xyMask);

      final BaseHash baseHash = BaseHashes.get(d0h);
      edgePixelNeighbours(d0h, this.fc.ij2i(hash), this.fc.ij2j(hash), neihbours);
      if (sorted) {
        neihbours.sortByHashAsc();
      }
      for(int i = 0; i < neihbours.size(); i++) {
        final long neigHash = neihbours.get(i);
        final MainWind neigDirection = neihbours.getDirection(i);
        appendSortedInternalEdgeElement(neigHash, toEdgeDeltaDepth,
            baseHash.getDirectionFromNeighbour(neigDirection), result);
      }
    } else {
      // Easy: always use opposite direction
      innerPixelNeighbours(d0hBits, iInD0hBits, jInD0hBits, neihbours);
      if (sorted) {
        neihbours.sortByHashAsc();
      }
      for(int i = 0; i < neihbours.size(); i++) {
        final long neigHash = neihbours.get(i);
        final MainWind neigDirection = neihbours.getDirection(i);
        appendSortedInternalEdgeElement(neigHash, toEdgeDeltaDepth,
            neigDirection.getOppositeDirection(), result);
      }
    }
  }

  private final void checkHashRange(long hash) {
    if (hash < 0 || this.nHash <= hash) {
      throw new IllegalArgumentException("Hash value " + hash + " must be in [0, " + this.nHash + "[");
    }
  }

  private boolean isInBasePixelBorderFromBits(
      final long iInBasePixelBits, final long jInBasePixelBits) {
    return 0 == iInBasePixelBits || iInBasePixelBits == this.xMask
        || 0 == jInBasePixelBits || jInBasePixelBits == this.yMask;
  }


  private void appendSortedInternalEdgeElement(final long hash, final int toEdgeDeltaDepth,
      final MainWind direction, final FlatHashList result) {
    if (direction.isCardinal()) {
      result.put(internalCorner(hash, toEdgeDeltaDepth, direction.toCardinal()));
    } else if (direction.isOrdinal()) {
      // We could have avoided array copies here!!
      result.put(sortedInternalEdge(hash, toEdgeDeltaDepth, direction.toOrdinal()));
    } else {
      throw new IllegalArgumentException("Main wind " + direction + " is neither ordinal not cradinal.");
    }
  }

  /*private void innerPixelNeighbours(long d0hBits, long iBits, long jBits, final FlatHashList result) {
    // Could have simply been:
    //     innerPixelNeighbours(d0hBits, iBits, jBits, EnumSet.allOf(MainWind.class), result);
    // but we preferred to unroll the for loop.
    long ij = this.fc.hash2ij(iBits | jBits);
    final int i = this.fc.ij2i(ij);
    final int j = this.fc.ij2j(ij);
    // Compute i-1 and j-1 bits.
    // Depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could result
    // in making fewer operations
    ij = this.fc.ij2hash(i - 1, j - 1);
    final long jm1Bits = ij & this.yMask;
    final long im1Bits = ij & this.xMask;
    // Compute i+1 and j+1 bits.
    // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
    // result in making fewer operations
    ij = this.fc.ij2hash(i + 1, j + 1);
    final long jp1Bits = ij & this.yMask;
    final long ip1Bits = ij & this.xMask;
    result.put(bits2hash(d0hBits, im1Bits, jm1Bits));
    result.put(bits2hash(d0hBits,   iBits, jm1Bits));
    result.put(bits2hash(d0hBits, ip1Bits, jm1Bits));
    result.put(bits2hash(d0hBits, im1Bits,   jBits));
    result.put(bits2hash(d0hBits, ip1Bits,   jBits));
    result.put(bits2hash(d0hBits, im1Bits, jp1Bits));
    result.put(bits2hash(d0hBits,   iBits, jp1Bits));
    result.put(bits2hash(d0hBits, ip1Bits, jp1Bits));
  }*/

  private void innerPixelNeighbours(long d0hBits, long iBits, long jBits, final NeighbourList result) {
    // Could have simply been:
    //     innerPixelNeighbours(d0hBits, iBits, jBits, EnumSet.allOf(MainWind.class), result);
    // but we preferred to unroll the for loop.
    long ij = this.fc.hash2ij(iBits | jBits);
    final int i = this.fc.ij2i(ij);
    final int j = this.fc.ij2j(ij);
    // Compute i-1 and j-1 bits.
    // Depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could result
    // in making fewer operations
    ij = this.fc.ij2hash(i - 1, j - 1);
    final long jm1Bits = ij & this.yMask;
    final long im1Bits = ij & this.xMask;
    // Compute i+1 and j+1 bits.
    // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
    // result in making fewer operations
    ij = this.fc.ij2hash(i + 1, j + 1);
    final long jp1Bits = ij & this.yMask;
    final long ip1Bits = ij & this.xMask;
    result.put(bits2hash(d0hBits, im1Bits, jm1Bits), MainWind.S);
    result.put(bits2hash(d0hBits,   iBits, jm1Bits), MainWind.SE);
    result.put(bits2hash(d0hBits, ip1Bits, jm1Bits), MainWind.E);
    result.put(bits2hash(d0hBits, im1Bits,   jBits), MainWind.SW);
    result.put(bits2hash(d0hBits, ip1Bits,   jBits), MainWind.NE);
    result.put(bits2hash(d0hBits, im1Bits, jp1Bits), MainWind.W);
    result.put(bits2hash(d0hBits,   iBits, jp1Bits), MainWind.NW);
    result.put(bits2hash(d0hBits, ip1Bits, jp1Bits), MainWind.N);
  }
  private void innerPixelNeighbours(long d0hBits, long iBits, long jBits, EnumSet<MainWind> directions,
      final NeighbourList result) {
    // Part of this code redundant with "innerPixelNeighbours"
    // to avoid creation of object containing 4 longs
    long ij = this.fc.hash2ij(iBits | jBits);
    final int i = this.fc.ij2i(ij);
    final int j = this.fc.ij2j(ij);
    // Compute i-1 and j-1 bits.
    // Depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could result
    // in making fewer operations
    ij = this.fc.ij2hash(i - 1, j - 1);
    final long jm1Bits = ij & this.yMask;
    final long im1Bits = ij & this.xMask;
    // Compute i+1 and j+1 bits.
    // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
    // result in making fewer operations
    ij = this.fc.ij2hash(i + 1, j + 1);
    final long jp1Bits = ij & this.yMask;
    final long ip1Bits = ij & this.xMask;
    for (final MainWind mw : directions) {
      final long nHash = bits2hash(d0hBits,
          mw.pickRightSouthToEastLongValue(im1Bits, iBits, ip1Bits),
          mw.pickRightSouthToWestLongValue(jm1Bits, jBits, jp1Bits));
      addNeighbour(nHash, MainWind.SE, result);
    }
  }

  /*private void edgePixelNeighbours(long hash, final FlatHashList result) {
    final int d0h = (int) (hash >> this.twiceDepth);
    hash = this.fc.hash2ij(hash & this.xyMask);
    edgePixelNeighbours(d0h, this.fc.ij2i(hash), this.fc.ij2j(hash), result);
  }*/

  private void edgePixelNeighbours(long hash, final NeighbourList result) {
    final int d0h = (int) (hash >> this.twiceDepth);
    hash = this.fc.hash2ij(hash & this.xyMask);
    edgePixelNeighbours(d0h, this.fc.ij2i(hash), this.fc.ij2j(hash), result);
  }

  /*private void edgePixelNeighbours(int baseCellHash, int iInBaseCell, int jInBaseCell,
      final FlatHashList result) {
    // Could have simply been edgePixelNeighbours(hash, EnumSet.allOf(MainWind.class) result)
    // but we prefered to unroll the for loop.
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.S ), result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.SE), result);
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.E ), result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.SW), result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.NE), result);
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.W ), result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.NW), result);
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.N ), result);
  }*/

  private void edgePixelNeighbours(int baseCellHash, int iInBaseCell, int jInBaseCell,
      final NeighbourList result) {
    // Could have simply been edgePixelNeighbours(hash, EnumSet.allOf(MainWind.class) result)
    // but we prefered to unroll the for loop.
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.S ), MainWind.S , result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.SE), MainWind.SE, result);
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.E ), MainWind.E , result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.SW), MainWind.SW, result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.NE), MainWind.NE, result);
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.W ), MainWind.W , result);
    addNeighbour        (neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.NW), MainWind.NW, result);
    addNeighbourIfExists(neighbour(baseCellHash, iInBaseCell, jInBaseCell, MainWind.N ), MainWind.N , result);
  }

  private void edgePixelNeighbours(long hash, final EnumSet<MainWind> directions,
      final NeighbourList neighbours) {
    final int d0h = (int) (hash >> this.twiceDepth);
    hash = this.fc.hash2ij(hash & this.xyMask);
    final int iInD0h = this.fc.ij2i(hash);
    final int jInD0h = this.fc.ij2j(hash);
    for (final MainWind direction : directions) {
      final long neigHash = neighbour(d0h, iInD0h, jInD0h, direction);
      addNeighbourIfExists(neigHash, direction, neighbours);
    }
  }

  private long neighbour(int baseCellHash, int iInBaseCell, int jInBaseCell, MainWind direction) {
    final BaseHash baseHash = BaseHashes.get(baseCellHash);
    iInBaseCell += direction.getOffsetSE();
    jInBaseCell += direction.getOffsetSW();
    final MainWind neigBaseCellDirection = getNeighbourBaseCellDirection(iInBaseCell, jInBaseCell);
    baseCellHash = iInBaseCell;
    iInBaseCell = baseHash.pickRightIndexOnNeighbourSouthToEastAxis(neigBaseCellDirection,
        iInBaseCell, jInBaseCell, this.nsideRemainderMask);
    jInBaseCell = baseHash.pickRightIndexOnNeighbourSouthToWestAxis(neigBaseCellDirection,
        baseCellHash, jInBaseCell, this.nsideRemainderMask);
    baseCellHash = baseHash.getNeighbour(neigBaseCellDirection);
    return (((long) baseHash.getValue()) << this.twiceDepth) | this.fc.ij2hash(iInBaseCell, jInBaseCell);
  }

  private MainWind getNeighbourBaseCellDirection(final int iInBaseCell, final int jInBaseCell) {
    final int neigBaseCellOffsetSE = neighbourBaseCellOffset(iInBaseCell);
    final int neigBaseCellOffsetSW = neighbourBaseCellOffset(jInBaseCell);
    return MainWind.getFromOffset(neigBaseCellOffsetSE, neigBaseCellOffsetSW);
  }

  private int neighbourBaseCellOffset(final int cooInBaseCell) {
    return (cooInBaseCell >> 31) | (cooInBaseCell >> this.depth);
  }


  private static void addNeighbour(final long neighbourHash, final FlatHashList result) {
    result.put(neighbourHash);
  }
  /*private static void addNeighbourIfExists(final long neighbourHash, final FlatHashList result) {
    if (neighbourHash >= 0) {
      addNeighbour(neighbourHash, result);
    }
  }*/
  private static void addNeighbour(final long neighbourHash, final MainWind direction,
      final NeighbourList neighbours) {
    neighbours.put(neighbourHash, direction);
  }
  private static void addNeighbourIfExists(final long neighbourHash, final MainWind direction,
      final NeighbourList neighbours) { // Branch version
    if (neighbourHash >= 0) {
      addNeighbour(neighbourHash, direction, neighbours);
    }
  }


  // UTILIY METHODS

  private static long xMask(final int depth) {
    return 0x5555555555555555L >>> (64 - (depth << 1));
  }

  private static long yMask(final int depth) {
    return 0xAAAAAAAAAAAAAAAAL >>> (64 - (depth << 1));
  }

  private static long xyMask(final int depth) {
    // return (1L << (depth << 1)) - 1;
    return (-1L) >>> (64 - (depth << 1));
  }

}
