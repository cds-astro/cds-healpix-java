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

import java.util.EnumMap;
import java.util.EnumSet;

import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.CompassPoint.MainWind;
import cds.healpix.CompassPoint.Ordinal;
import cds.healpix.fillingcurve.FillingCurve2D;
import cds.healpix.fillingcurve.FillingCurve2DType;

import static cds.healpix.Healpix.nside;
import static cds.healpix.Healpix.getBestStartingDepth;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.toBits;
import static cds.healpix.common.math.Math.FOUR_OVER_PI;

/**
 * Implementation of the HEALPix tessellation using the NESTED scheme: the NESTED scheme consists
 * in the concatenation of 12 z-order curves, one by HEALPix base cell.<br>
 * The hash value in the NESTED scheme is build as follow:
 * <b> 0...0bbbb112233...</b><br>
 * With<br>
 * <ul>
 *   <li>0...0: unused bits</li>
 *   <li>bbbb: the 4 bits coding the base cell, in [0, 11]</li>
 *   <li>11: the 2 bits of the z-order curve coding depth 1</li>
 *   <li>22: the 2 bits of the z-order curve coding depth 2</li> 
 *   <li>33: the 2 bits of the z-order curve coding depth 3</li> 
 *   <li>...</li> 
 * </ul>
 * 
 * 
 * The constructor is {@code protected}. to get an instance, see the {@link Healpix#getNested(int)}
 * or {@link Healpix#getNested(int, FillingCurve2DType)}.
 *
 * @author F.-X. Pineau
 *
 */
public final class HealpixNested implements HashComputer, VerticesAndPathComputer, NeighbourSelector {

  /** Lookup table to find the base-cell identifier from the rotated and shifted projection.
   *  Tests show better performances with respect to:
   *   j = 4 - (j + i);
   *   i = ((i - (j >>> 63)) & 3) + (++j << 2)
   * But it may depends on the hardware and on the cache occupation.
   */
  /*static final byte[][] D0C_LOOKUP = new byte[][] {
    {-1, -1, -1,  8,  4}, //   ----> y-axis 
    {-1, -1,  9,  5,  0}, //  |
    {-1, 10,  6,  1},     //  |
    {11,  7,  2},         //  v
    { 4,  3,},            // x-axis
  };*/
  
  /** Lookup tables to retrieve the bits of the base (depth=0) cell at all possible depth. */
  // private static final long[][][] D0C_BITS_LUPTS = new long[30][][];
  
  /** Lookup table to retrieve the bits of the base (depth=0) cell at the object depth .
   *  We modified {@link #D0C_LOOKUP} to take into account cases of lat = +-90 deg 
   *  and equatorial/polar caps limit.
   * */
  // final long[][] d0cBitsLUPT;
  
  final int depth;

  final long nHash;
  
  // Derived quantities to speed up computations
  final int twiceDepth;
  final int nside;
  final double oneOverNside;
  private final long halfNside4IEEEdouble;
  
  final int nsideRemainderMask;
  private final double xScale;    // [0,  pi/4] --> [0, nside]

  private final double yScaleCEA; // [0, 2 / 3] --> [0, nside]
  private final double yScaleCOL; // [1,     0] --> [nside, 2*nside]
  final double xHalfScale;
  final double yHalfScaleCEA;
  final double yHalfScaleCOL;

  /** Mask used to retrieve the value of the base (i.e., depth=0) pixel in
   * the full hash. So all bits are : 15 << (this.depth << 1).
   * Background: we code 2^n - 1 values on n bits */
  private final long d0Mask;
  /** Mask used to retrieve the yx part of a hash, i.e. removing the value
   * of the base pixel: the 2*nside least significant bits are set to 1.  */
  final long xyMask;
  /** Mask used to retrieve the bits of the x value only in the full hash. */
  final long xMask;
  /** Mask used to retrieve the bits of the y value only in the full hash. */
  final long yMask;
  
  final FillingCurve2DType fillingCurveType;
  final FillingCurve2D fc;
  // private final FillingCurve2D fc2Unshuflle;
  
  private HashComputer hashComputer;
  private VerticesAndPathComputer vpComputer;
  private HealpixNestedNeighbourSelector neigSelector;

  protected HealpixNested(final int depth, final FillingCurve2DType fillingCurveType) {
    this.depth = depth;
    this.nside = nside(this.depth);
    // this.d0cBitsLUPT = getD0cBitsLookupTable(this.depth);
    this.halfNside4IEEEdouble = Healpix.halfNside4IEEEdouble(this.depth);
    this.oneOverNside = 1 / (double) this.nside;
    this.nHash = Healpix.nHash(this.depth);
    this.twiceDepth = this.depth << 1;           assert twiceDepth == this.depth * 2;
    this.d0Mask = 15L << this.twiceDepth; // 15: 0...01111 in binary
    if (this.depth > 0) {
      this.xyMask = (1L << this.twiceDepth) - 1;
      this.xMask = 0x5555555555555555L >>> (64 - this.twiceDepth); // ...0101
      this.yMask = this.xMask << 1;
      assert this.yMask == 0xAAAAAAAAAAAAAAAAL >>> (64 - this.twiceDepth); // ...1010
    } else {
      this.xyMask = 0;
      this.xMask = 0;
      this.yMask = 0;
    }
    this.nsideRemainderMask = this.nside - 1;
    this.xScale = this.nside * FOUR_OVER_PI;
    this.yScaleCEA = this.nside * 1.5d; // i.e. * 3 / 2
    this.yScaleCOL = this.nside;
    this.xHalfScale = 0.5 * this.xScale;
    this.yHalfScaleCEA = 0.5 * this.yScaleCEA;
    this.yHalfScaleCOL = 0.5 * this.yScaleCOL;
    this.fillingCurveType = fillingCurveType;
    this.fc = this.fillingCurveType.get(depth);
    this.neigSelector = new HealpixNestedNeighbourSelector(this);
  }

  /*private static long[][] getD0cBitsLookupTable(final int depth) {
    long[][] ref = D0C_BITS_LUPTS[depth];
    if (ref == null) {
      synchronized(D0C_BITS_LUPTS) {
        ref = D0C_BITS_LUPTS[depth];
        if (ref == null) {
          ref = buildLookUpTable(depth);
          D0C_BITS_LUPTS[depth] = ref;
        }
      }
    }
    return ref;
  }*/
  
  /*private static long[][] buildLookUpTable(int depth) {
    final int twiceDepth = depth << 1;
    long xMask = 0;
    long yMask = 0;
    long xyMask = 0;
    if (depth > 0) {
      xyMask = (1L << twiceDepth) - 1;
      xMask = 0x5555555555555555L >>> (64 - twiceDepth); // ...0101
      yMask = xMask << 1;                                // ...1010
    }
    final long nill = -1L;
    final long bc00 =  0L << twiceDepth;
    final long bc01 =  1L << twiceDepth;
    final long bc02 =  2L << twiceDepth;
    final long bc03 =  3L << twiceDepth;
    final long bc04 =  4L << twiceDepth;
    final long bc05 =  5L << twiceDepth;
    final long bc06 =  6L << twiceDepth;
    final long bc07 =  7L << twiceDepth;
    final long bc08 =  8L << twiceDepth;
    final long bc09 =  9L << twiceDepth;
    final long bc10 = 10L << twiceDepth;
    final long bc11 = 11L << twiceDepth;
    final long mg0y = (0L << twiceDepth) | yMask;
    final long mg1y = (1L << twiceDepth) | yMask;
    final long mg2y = (2L << twiceDepth) | yMask;
    final long mg3y = (3L << twiceDepth) | yMask;
    final long m0xy = (0L << twiceDepth) | xyMask;
    final long m1xy = (1L << twiceDepth) | xyMask;
    final long m2xy = (2L << twiceDepth) | xyMask;
    final long m3xy = (3L << twiceDepth) | xyMask;
    return new long[][] {
      {nill, nill, nill, bc08, bc04},       //   ----> y-axis 
      {nill, nill, bc09, bc05, bc00, mg0y}, //  |
      {nill, bc10, bc06, bc01, mg1y, m0xy}, //  |
      {bc11, bc07, bc02, mg2y, m1xy},       //  v
      {bc04, bc03, mg3y, m2xy},             // x-axis
      {nill, mg0y, m3xy}
    };
  }*/
  
  
  
  /*
   * Project the given position on the unit sphere onto an Euclidean two-dimensional plane using
   * the HEALPix projection of center (0, 0) and parameters H=4 anf K=3 (see Calabrett2007).
   * The scale of the projection with respect to the canonical projection is 4*nside/pi so that
   * the coordinates of each pixel center are integers.
   * @param lonlat position on the unit sphere to be projected.
   * @return the projected coordinates, X in [-4*nside, 4*nside] and Y in [-2*nside, 2*nside],
   *         so that each pixel center and vertices have integer coordinates. 
   */
  /*public XY project(final LonLat lonlat) {
    final SettableXY result = new SettableXYImpl();
    this.project(lonlat, result);
    return result;
  }*/

  /*
   * Same function as {@link #proj(LonLat)} but avoiding the instantiation of the {@link XY}
   * object. One possible use is for example to create one {@link SetableXY} object by thread and to
   * compute successive projections (in a given thread) by reusing the {@link SetableXY} object.
   * @param lonlat position on the unit sphere to be projected.
   * @param result projected coordinates onto the Euclidean plane.
   */
 /* public void project(final LonLat lonlat, final SettableXY result) {
    if (Healpix.isInEquatorialRegion(lonlat)) {
      this.projectInCylindricalEquaArea(lonlat, result);
    } else {
      this.projectInCollignon(lonlat, result);
    }
  }

  private void projectInCylindricalEquaArea(final LonLat lonlat, final SettableXY result) {
    Healpix.projCylindricalEquaAreaNoScale(lonlat, result);
    result.setX(result.x() * this.xScale);
    result.setY(result.y() * this.yScaleCEA);
  }

  private void projectInCollignon(final LonLat lonlat, final SettableXY result) {
    Healpix.projCollignonNoScale(lonlat, result);
    result.setX(result.x() * this.xScale);
    result.setY(result.y() * this.yScaleCOL);
  }*/


  /*
   * Perform the inverse operation of the projection, {@see #proj(LonLat)}.
   * @param xy x must be in [-4*nside, 4*nside] and y must be in [-2*nside, 2*nside]
   * @return x must be in [-4*nside, 4*nside] and y must be in [-2*nside, 2*nside]
   */
  /*public LonLat deProject(final XY xy) {
    return Healpix.deproj(xy.x() * this.scaleInv, xy.y() * this.scaleInv);
  }*/

  /*
   * Perform the inverse operation of the projection, {@see #proj(LonLat)}.
   * @param xy x must be in [-4*nside, 4*nside] and y must be in [-2*nside, 2*nside]
   * @param result the result of the deprojection
   */
  /*public void deProject(final XY xy, final SettableLonLat result) {
    Healpix.deproj(xy.x() * this.scaleInv, xy.y() * this.scaleInv, result);
  }*/


  @Override 
  public int depth() {
    return this.depth;
  }
  
  /**
   * In multi-threaded environments, each thread must have its own {@link HashComputer} since
   * {@link HashComputer} is not thread safe.
   * @return a new {@link HashComputer}, to be used in a distinct thread.
   */
  public HashComputer newHashComputer() {
    return new HealpixNestedHashComputer(this);
  }

  /**
   * WARNING: the return object in not thread-safe!
   * @param auxAxis object describing the auxiliary axis
   * @return a new {@link HashComputerWithAux}, to be used in a distinct thread.
   */
  public HashComputerWithAux newHashComputerWithAux(final AuxiliaryAxis auxAxis) {
    return new HealpixNestedHashComputerWithAux(this, auxAxis);
  }
  
  
  /**
   * This method is thread safe. For better performances, use one {@link HashComputer} created by
   * {@link #newHashComputer()} per thread.
   * WARNING: this method IS NOT thread safe. In multi-threaded environments, use the method
   * {@link #newHashComputer()} to create a {@link HashComputer} for each thread.
   * @param lonRad longitude (in radians) of the the coordinate we want the HEALPix hash value.
   * @param latRad latitude (in radians) of the coordinate we want the HEALPix hash value.
   * @return the HEALPix hash value of the given coordinate.
   */
  @Override
  public long hash(final double lonRad, final double latRad) {
    if (this.hashComputer == null) {
      this.hashComputer = new HealpixNestedHashComputer(this);
    }
    return this.hashComputer.hash(lonRad, latRad);
  }

  /**
   * Returns a {@link VerticesAndPathComputer}, to be used in a distinct thread.
   * @return a {@link VerticesAndPathComputer}, to be used in a distinct thread.
   */
  public VerticesAndPathComputer newVerticesAndPathComputer() {
    return new HealpixNestedVerticesAndPathComputer(this);
  }
  
  @Override
  public double[] center(long hash) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    return this.vpComputer.center(hash);
  }

  @Override
  public void center(long hash, double[] resultLonLat) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    this.vpComputer.center(hash, resultLonLat);
  }

  @Override
  public double[] vertex(long hash, Cardinal cardinalPoint) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    return this.vpComputer.vertex(hash, cardinalPoint);
  }

  @Override
  public void vertex(long hash, Cardinal cardinalPoint, double[] resultLonLat) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    this.vpComputer.vertex(hash, cardinalPoint, resultLonLat);
  }

  @Override
  public EnumMap<Cardinal, double[]> vertices(long hash, EnumSet<Cardinal> cardinalPoints) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    return this.vpComputer.vertices(hash, cardinalPoints);
  }

  @Override
  public void vertices(long hash, EnumMap<Cardinal, double[]> cardinalPoints) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    this.vpComputer.vertices(hash, cardinalPoints);
  }

  @Override
  public double[][] pathAlongCellSide(long hash, Cardinal fromVertex, Cardinal toVertex,
      boolean isToVertexIncluded, int nSegments) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    return this.vpComputer.pathAlongCellSide(hash, fromVertex, toVertex, isToVertexIncluded, nSegments);
  }

  @Override
  public void pathAlongCellSide(long hash, Cardinal fromVertex, Cardinal toVertex,
      boolean isToVertexIncluded, int nSegments, double[][] pathPoints) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    this.vpComputer.pathAlongCellSide(hash, fromVertex, toVertex, isToVertexIncluded, nSegments, pathPoints);
  }

  @Override
  public double[][] pathAlongCellEdge(long hash, Cardinal startingVertex,
      boolean clockwiseDirection, int nSegmentsBySide) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    return this.vpComputer.pathAlongCellEdge(hash, startingVertex, clockwiseDirection, nSegmentsBySide);
  }

  @Override
  public void pathAlongCellEdge(long hash, Cardinal startingVertex, boolean clockwiseDirection,
      int nSegmentsBySide, double[][] pathPoints) {
    if (this.vpComputer == null) {
      this.vpComputer = new HealpixNestedVerticesAndPathComputer(this);
    }
    this.vpComputer.pathAlongCellEdge(hash, startingVertex, clockwiseDirection, nSegmentsBySide, pathPoints);
  }
  
  /**
   * The internal {@link NeighbourSelector} being thread-safe and unmutable, we return it instead
   * of creating a new instance.
   * @return the internal neighbour selector.
   */
  public NeighbourSelector newNeighbourSelector() {
    return this.neigSelector;
  }
  
  @Override
  public long neighbour(long hash, MainWind direction) {
    return this.neigSelector.neighbour(hash, direction);
  }

  @Override
  public NeighbourList neighbours(long hash) {
    return this.neigSelector.neighbours(hash);
  }

  @Override
  public void neighbours(long hash, NeighbourList result) {
    this.neigSelector.neighbours(hash, result);
  }

  @Override
  public void neighbours(long hash, FlatHashList result) {
    this.neigSelector.neighbours(hash, result);
  }

  @Override
  public NeighbourList neighbours(long hash, EnumSet<MainWind> directions) {
    return this.neigSelector.neighbours(hash, directions);
  }

  @Override
  public void neighbours(long hash, EnumSet<MainWind> directions, NeighbourList result) {
    this.neigSelector.neighbours(hash, directions, result);
  }

  @Override
  public FlatHashList internalEdges(long hash, int deltaDepth) {
    return this.neigSelector.internalEdges(hash, deltaDepth);
  }

  @Override
  public void internalEdges(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.internalEdges(hash, deltaDepth, result);
  }

  @Override
  public FlatHashList sortedInternalEdges(long hash, int deltaDepth) {
    return this.neigSelector.sortedInternalEdges(hash, deltaDepth);
  }

  @Override
  public void sortedInternalEdges(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.sortedInternalEdges(hash, deltaDepth, result);
  }

  @Override
  public FlatHashList sortedInternalEdge(long hash, int deltaDepth, Ordinal direction) {
    return this.neigSelector.sortedInternalEdge(hash, deltaDepth, direction);
  }

  @Override
  public void sortedInternalEdge(long hash, int deltaDepth, Ordinal direction,  FlatHashList result) {
    this.neigSelector.sortedInternalEdge(hash, deltaDepth, direction, result);
  }

  @Override
  public FlatHashList sortedInternalEdgeNE(long hash, int deltaDepth) {
    return this.neigSelector.sortedInternalEdgeNE(hash, deltaDepth);
  }

  @Override
  public void sortedInternalEdgeNE(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.sortedInternalEdgeNE(hash, deltaDepth, result);
  }

  @Override
  public FlatHashList sortedInternalEdgeNW(long hash, int deltaDepth) {
    return this.neigSelector.sortedInternalEdgeNW(hash, deltaDepth);
  }

  @Override
  public void sortedInternalEdgeNW(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.sortedInternalEdgeNW(hash, deltaDepth, result);
  }

  @Override
  public FlatHashList sortedInternalEdgeSE(long hash, int deltaDepth) {
    return this.neigSelector.sortedInternalEdgeSE(hash, deltaDepth);
  }

  @Override
  public void sortedInternalEdgeSE(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.sortedInternalEdgeSE(hash, deltaDepth, result);
  }

  @Override
  public FlatHashList sortedInternalEdgeSW(long hash, int deltaDepth) {
    return this.neigSelector.sortedInternalEdgeSW(hash, deltaDepth);
  }

  @Override
  public void sortedInternalEdgeSW(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.sortedInternalEdgeSW(hash, deltaDepth, result);
  }

  @Override
  public long internalCorner(long hash, int deltaDepth, Cardinal direction) {
    return this.neigSelector.internalCorner(hash, deltaDepth, direction);
  }

  @Override
  public long internalCornerN(long hash, int deltaDepth) {
    return this.neigSelector.internalCornerN(hash, deltaDepth);
  }

  @Override
  public long internalCornerS(long hash, int deltaDepth) {
    return this.neigSelector.internalCornerS(hash, deltaDepth);
  }

  @Override
  public long internalCornerE(long hash, int deltaDepth) {
    return this.neigSelector.internalCornerE(hash, deltaDepth);
  }

  @Override
  public long internalCornerW(long hash, int deltaDepth) {
    return this.neigSelector.internalCornerW(hash, deltaDepth);
  }

  @Override
  public FlatHashList externalEdges(long hash, int deltaDepth) {
    return this.neigSelector.externalEdges(hash, deltaDepth);
  }

  @Override
  public void externalEdges(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.externalEdges(hash, deltaDepth, result);
  }

  @Override
  public FlatHashList sortedExternalEdges(long hash, int deltaDepth) {
    return this.neigSelector.sortedExternalEdges(hash, deltaDepth);
  }

  @Override
  public void sortedExternalEdges(long hash, int deltaDepth, FlatHashList result) {
    this.neigSelector.sortedExternalEdges(hash, deltaDepth, result);
  }
  
  /**
   * Returns a {@link HealpixNestedFixedRadiusConeComputer}, to be used in a distinct thread.
   * @param coneRadiusRad radius of the cone.
   * @return a {@link HealpixNestedFixedRadiusConeComputer}, to be used in a distinct thread.
   */
  public HealpixNestedFixedRadiusConeComputer newConeComputer(final double coneRadiusRad) {
    if (coneRadiusRad >= Math.PI) {
      return new NestedAllSky(coneRadiusRad, this.depth);
    }
    final int optimalStartingDepth = getBestStartingDepth(coneRadiusRad);
    if (optimalStartingDepth == -1) { // r in [48deg, PI[
      return new NestedSmallCellApproxedMethod(optimalStartingDepth, this.depth, coneRadiusRad);
    } else if (this.depth <= optimalStartingDepth) {
      if (coneRadiusRad < 1e-6) { // ~ 200 mas
        return new NestedLargeCellSmallRadius(this, coneRadiusRad);
      } else {
        return new NestedLargeCell(this, coneRadiusRad);
      }
    } else if (coneRadiusRad < 1e-6) { // ~ 200 mas
      return new NestedSmallCell(optimalStartingDepth, this.depth, coneRadiusRad, SmallConeOrdinalHashComputer.UI);
    } else {
      return new NestedSmallCell(optimalStartingDepth, this.depth, coneRadiusRad, RegularConeOrdinalHashComputer.UI);
    }
  };

  /**
   * Simpler version of the cone search, but at the cost of approximation (false positieves).
   * @param coneRadiusRad radius of the cone
   * @return an object allowing to perform fixed radius cone search queries.
   */
  public HealpixNestedFixedRadiusConeComputer newConeComputerApprox(final double coneRadiusRad) {
    if (coneRadiusRad >= Math.PI) {
      return new NestedAllSky(coneRadiusRad, this.depth);
    }
    final int optimalStartingDepth = getBestStartingDepth(coneRadiusRad);
    if (optimalStartingDepth >= this.depth) {
      return new NestedLargeCellApproxedMethod(optimalStartingDepth, this.depth, coneRadiusRad);
    } else {
      return new NestedSmallCellApproxedMethod(optimalStartingDepth, this.depth, coneRadiusRad);
    }
  }
  
  /**
   * We suppose here that the radius of the (cross-match) cone is smaller that the typical size of
   * the cell by a factor at least 2 or 3.
   * The algorithm is design to be very fast, at the cost of very raw approxiations (false positives).
   * @param coneRadiusRad radius of the cone
   * @return an object allowing to perform fixed radius cone search queries.
   */
  public HealpixNestedFixedRadiusCone4XMatch newConeComputer4Xmatch(final double coneRadiusRad) {
    return new HealpixNestedFixedRadiusCone4XMatch(this.depth, coneRadiusRad);
  }
  
  /**
   * Returns a new instance of {@link HealpixNestedPolygonComputer} at this object depth.
   * @return a new instance of {@link HealpixNestedPolygonComputer} at this object depth.
   */
  public HealpixNestedPolygonComputer newPolygonComputer() {
    return new HealpixNestedPolygonComputer(this);
  }
  
  /**
   * Returns the ring representation if the given nested hash.
   * @param nestedHash an HEALPix cell nested hash.
   * @return the ring representation if the given nested hash.
   */
  public long toRing(final long nestedHash) {
   final HashParts hparts = decodeRegularRing(nestedHash);
   return regularRingEncode(hparts.baseCellHash(), hparts.iInBaseCell(), hparts.jInBaseCell());
  }

  /**
   * Encode the given data into the ring scheme hash.
   * @param baseCellHash base cell number, in [0, 12[.
   * @param iInBasePixel x-axis (i.e. South-East) z-curve coordinate in the base cell.
   * @param jInBasePixel y-axis (i.e. South West) z-curve coordinate in the base cell.
   * @return the ring encoded hash of the given base cell and coordinates in the base cell.
   */
  public long regularRingEncode(long baseCellHash, int iInBasePixel, int jInBasePixel) {
    /*final int depth = 2;
    final HealpixNested hn = Healpix.getNested(depth);
    final List<String> lines = new ArrayList<String>((int) hn.nHash + 1);
    lines.add("nested,ring,ra,dec,h,l,iRing,firstIsoLatIndex,iInRing");
    for (int i = 0; i < hn.nHash; i++) {
      final HashParts hparts = hn.decodeRegularHash(i);
      final double[] center = hn.newVerticesAndPathComputer().center(i);  
      baseCellHash = hparts.baseCellHash();
      iInBasePixel = hparts.iInBaseCell();
      jInBasePixel = hparts.jInBaseCell();*/
      
      
    // Number of isolatitude rings in a cell = hmax = 2*nside - 1 => index max = 2*(nside - 2)
    // => index of ring at lat=0 = 2*(2*nside - 1) = 2*nside - 1 => always odd
    // h = x + y (<=> rotation 45 and scale sqrt(2))
    // North polar cap:   i =  hmax - h            = 2*nside - 2 - (x+y)
    // Equatorial region: i = (hmax - h) +   nside = 3*nside - 2 - (x+y)
    // South polar cap:   i = (hmax - h) + 2*nside = 4*nside - 2 - (x+y)
    final long h = iInBasePixel + jInBasePixel;
    final long l = iInBasePixel - jInBasePixel;
    final long jBasePixel = dividedBy4Quotient(baseCellHash);    assert 0 <= jBasePixel && jBasePixel <= 2;
    final long iRing = nsideTime(jBasePixel + 2) - (h + 2);      assert 0 <= iRing  &&  iRing <  Healpix.nIsolatitudeRings(this.depth);
    // Number of elements in isolatitude ring of index i (i in [0, 2*(2*nside - 1)]):
    // North polar cap: if (i < nside)         nj = 4 * i
    // South polar cap: if (i >= 3*nside - 1)  nj = 4 * h = 4*((4*nside-2) - i) 
    // Equatorial regi: if (ns <= i < 3*ns-1)  nj = 4 * nside
    // l = x - y; In a base cell, l in [-nside+1, nside-1] => 2*nside - 1 values
    // EQR: j = l / 2 + nside/2 (if not equatorial cell) + nside*(ipix%4) (special case if l<0 && baseCell=4)
    // NPC: j = l / 2 + (i + 1) / 2 + (i+1)*(ipix%4)
    // SPC: j = l / 2 + (h + 1) / 2 + (h+1)*(ipix%4)
    long firstIsolatIndex = 0;
    long iInRing = dividedBy2Quotient(l); // Such that: 1/2 = 0; -1/2 = -1
    if (iRing < nside) { // North polar cap
      // sum from i = 1 to ringIndex of 4 * i = 4 * i*(i+1)/2 = 2 * i*(i+1)
      final long ip1 = iRing + 1;
      firstIsolatIndex = (iRing * ip1) << 1;
      iInRing += dividedBy2Quotient(ip1);
      iInRing += ip1 * modulo4(baseCellHash);
    } else if (iRing >= nsideTime(3) - 1) { // South polar cap
      final long ip1 = h + 1;
      firstIsolatIndex = nHash - twiceValTimeValPlusOne(ip1); //((ip1 * (ip1 + 1)) << 1);
      iInRing += dividedBy2Quotient(ip1);
      iInRing += ip1 * modulo4(baseCellHash);
    } else { // Equatorial region: sum of i from i=1 to nside
      // firstIsolatIndex = ((nside  * (nside + 1)) << 1) + ((ringIndex - nside) * (nside << 2));
      firstIsolatIndex = twiceNsideTimeNsidePlusOne() + minusNsideTimeFourNside(iRing); // ((ringIndex - nside) * (nside << 2));
      iInRing += nsideTime(modulo2(jBasePixel + 1)) >> 1;
      iInRing += nsideTime(baseCellHash == 4 && l < 0 ? 4 : modulo4(baseCellHash));
    }
    return firstIsolatIndex + iInRing;
    
    // From the paper:
    // (v = vertical index = x+y)
    // f = base resolution number
    // fr = I(f/N_LON) with N_LON = 4
    // F1(f) = fr + 2
    // long pixelInRingIndex = ( F2(f) Nside + l + s) / 2;
    // l = x - y
    // F2(f) = 2 * (f mod 4) - (fr mod 2) + 1
    // s = phase shift = (i - Nside + 1) mod 2 in North equtorial belt and s=1 in North polar cap
    // EQR: s = (i-nside+1)%2; F2(f) = ipix%4    ; l / 2 + nside*(ipix%4) + (i-1)%2;
    // NPC: s = 1            ; F2(f) = ipix%4 + 1; l / 2 + nside*(ipix%4) + nside / 2 + 1 / 2
    // SPC: s = 1            ; F2(f) = ipix%4 + 1; l / 2 + nside*(ipix%4) + nside / 2 + 1 / 2
    
    
      // long ring = hn.regularRingEncode(hparts.baseCellHash(), hparts.iInBaseCell(), hparts.jInBaseCell());
      /*lines.add(i + "," + (firstIsolatIndex + iInRing) + "," + String.format(Locale.US, "%.6f,%.6f,%d,%d,%d,%d,%d",
          center[0], center[1], h, l, iRing, firstIsolatIndex, iInRing));

    }
    final File f = new File("./nested2ring.depth" + depth + ".csv");
    Files.write(f.toPath(), lines, StandardCharsets.US_ASCII);
    return 0;*/
  }

  /**
   * Returns the nested representation if the given ring hash.
   * @param ringHash an HEALPix cell ring hash.
   * @return the nested representation if the given ring hash.
   */
  public long toNested(long ringHash) {
    final HashParts hparts = decodeRegularRing(ringHash);
    return encodeRegularHash(hparts.baseCellHash(), hparts.iInBaseCell(), hparts.jInBaseCell());
  }
  
  HashParts decodeRegularRing(long hash) { // pullApart
    // 4*sum from i=1 to nside of i =  4 * nside(nside+1)/2 = 2*nside*(nside+1)
    final long firstHashInEquatorialRegion = twiceNsideTimeNsidePlusOne();
    final long firstHashInSouthPolarCap = nHash - firstHashInEquatorialRegion;
    long iRing = 0, iInRing = 0;
    int baseCellHash;
    long h, l, iInBaseCell, jInBaseCell;
    if (hash < firstHashInEquatorialRegion) { // North polar cap
      // 2*n(n+1) = x => 2n^2+2n-x = 0 => b^2-4ac = 4+8x = 4(1+2x)
      //              => n = [-2+2*sqrt(1+2x)]/4 => n = [sqrt(1+2x) - 1] / 2
      // => n^2+n-x/2 = 0 => b^2-4ac = 1 + 2x => n = [sqrt(1+2x) - 1] / 2
      // n - 1 = ring index
      // See here for a fast SQTR function: http://atoms.alife.co.uk/sqrt/index.html
      // but only for integer, not long :o/
      iRing = (((long) Math.sqrt(1 + (hash << 1))) - 1) >>> 1;
      final long nInRing = iRing + 1;
      iInRing = hash - twiceValTimeValPlusOne(iRing);
      baseCellHash = (int) (iInRing / nInRing);
      h = ((nside << 1) - 2) - iRing;
      l = ((iInRing - nInRing * baseCellHash) << 1) - iRing;
      iInBaseCell = (h + l) >>> 1;
      jInBaseCell = (h - l) >>> 1;
    } else if (hash >= firstHashInSouthPolarCap) { // South polar cap
      hash = nHash - 1 - hash; // start counting in reverse order from south polar cap
      iRing = (((long) Math.sqrt(1 + ((hash) << 1))) - 1) >> 1;
      final long nInRing = iRing + 1;
      // nElemsOnRing = ringIndex << 2;
      iInRing = ((nInRing << 2) - 1) - (hash - twiceValTimeValPlusOne(iRing));
      baseCellHash = (int) (iInRing / nInRing);
      h = iRing;
      l = ((iInRing - nInRing * baseCellHash) << 1) - iRing;
      baseCellHash += 8;
      iInBaseCell = (h + l) >>> 1;
      jInBaseCell = (h - l) >>> 1;
    } else { // Equatorial region
      // Set origin of ring indexes at the center of small cell in north corner of base cell 4 (Noth to South direction)
      iRing = hash - firstHashInEquatorialRegion; 
      iInRing = iRing;
      // <=> /= 4*nside (number of hash per line) => count the number of line from first equatorial line
      iRing >>= depth + 2;
      // Subtract number of hash in previous rings (-= nRings * 4*nside)
      iInRing -= iRing << (depth + 2);
      l = (iInRing << 1) + modulo2(iRing);
      // Set origin of h axis at center of small cell in south corner of base cell 4 (South to North direction)
      h = ((nside << 1) - 2) - iRing;
      // Rotation of -45
      iInBaseCell = (h + l) >>> 1;
      jInBaseCell = (h - l) >>> 1;
      // Offset of 4*nside in j
      jInBaseCell += nside << 2;
      int iBaseCell = dividedByNsideQuotient(iInBaseCell);
      int jBaseCell = dividedByNsideQuotient(jInBaseCell);
      baseCellHash = getBaseCellHash(iBaseCell, jBaseCell);
      iInBaseCell = moduloNside(iInBaseCell);
      jInBaseCell = moduloNside(jInBaseCell);
    }
    return new HashPartsImpl(baseCellHash, (int) iInBaseCell, (int) jInBaseCell);
  }
  
  private static int getBaseCellHash(int i, int j) {
    // retrun D0C_LOOKUP[iBaseCell][jBaseCell]; 
    // Tests on new computer show that NOBRANCH performs netter than LOOKUP_MATRIX
    // (and uses less CPU cache!).
    j = 4 - (j + i);
    return ((i - (j >>> 63)) & 3) + (++j << 2);
  }
  
  private long minusNsideTimeFourNside(long l) {
    // ((ringIndex - nside) * (nside << 2))
    return (l - nside) << (depth + 2);
  }
  
  private long twiceValTimeValPlusOne(long val) {
    return (val * (val + 1)) << 1;
  }
  
  private long twiceNsideTimeNsidePlusOne() {
    //   2*nside*(nside + 1)
    // = 2*[nside^2 + nside]
    return ((1L << (depth << 1)) + nside) << 1; 
  }

  
  double timeHalfNsideP(double v) {
    assert v >= 0;
    return fromBits(toBits(v) + this.halfNside4IEEEdouble);
  }
  /*private double divideByNsideP(double v) {
    assert v >= 0;
    return fromBits(toBits(v) - this.nside4IEEEdouble);
  }*/
  
  
  public int nsideTime(int i) {
    return i << depth;
  }

  public long nsideTime(long i) {
    return i << depth;
  }
  
  private long dividedBy2Quotient(final long d) {
    return d >> 1;
  }
  
  public long dividedBy4Quotient(final long d) {
    return d >> 2;
  }

  private int dividedByNsideQuotient(final int d) {
    return d >> this.depth;
  }

  int dividedByNsideQuotient(final long d) {
    return (int) (d >> this.depth);
  }

  private int dividedByNsideRemainder(final int d) {
    assert (d & this.nsideRemainderMask) == (d % this.nside);
    return d & this.nsideRemainderMask;
  }

  long modulo2(final long d) {
    return d & 1;
  }
  
  long modulo4(final long d) {
    return d & 3;
  }
  
  int moduloNside(final int d) {
    return d & this.nsideRemainderMask;
  }

  int moduloNside(final long d) {
    return ((int) d) & this.nsideRemainderMask;
  }

  static boolean isInEquatorialRegion(final int basePixel) {
    assert 0 <= basePixel &&  basePixel < 12;
    return (basePixel | 4) == 4;
  }

  static boolean isInNorhtPolarCapRegion(final int basePixel) {
    assert 0 <= basePixel &&  basePixel < 12;
    return (basePixel | 12) == 0;
  }

  static boolean isInSouthPolarCapRegion(final int basePixel) {
    assert 0 <= basePixel &&  basePixel < 12;
    return (basePixel | 8) == 8;
  }

  

  long encodeRegularHash(long basePixel, int iInBasePixel, int jInBasePixel) {
    return (basePixel << this.twiceDepth) | this.fc.ij2hash(iInBasePixel, jInBasePixel);
  }

  static long bits2hash(long basePixelBits, long iInBasePixelBits, long jInBasePixelBits) {
    return basePixelBits | iInBasePixelBits | jInBasePixelBits;
  }

  HashParts decodeRegularHash(long hash) { // pullApart
    final int d0h = (int) (hash >> this.twiceDepth);
    hash = this.fc.hash2ij(hash & this.xyMask);
    return new HashPartsImpl(d0h, this.fc.ij2i(hash), this.fc.ij2j(hash));
  }

  /**
   * Extract the value of the base pixel and the coordinates in the base pixel from the given hash.
   * @param hash
   * @param result
   */
  void decodeRegularHash(long hash, final SettableHashParts result) {
    result.setBaseCellHash((int) (hash >> this.twiceDepth));
    hash = this.fc.hash2ij(hash & this.xyMask);
    result.setIInBaseCell(this.fc.ij2i(hash));
    result.setJInBaseCell(this.fc.ij2j(hash));
  }

  HashBits pullBitsApart(long hash) {
    return new HashBits(hash & this.d0Mask, hash & this.xMask, hash & this.yMask);
  } 

}
