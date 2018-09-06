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


import java.util.EnumSet;

import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.CompassPoint.MainWind;
import cds.healpix.CompassPoint.Ordinal;
import cds.healpix.fillingcurve.FillingCurve2D;

import static cds.healpix.HealpixNested.bits2hash;
import static cds.healpix.common.math.Math.pow;

final class HealpixNestedNeighbourSelector implements NeighbourSelector {
  // should be an HealpixNested inner class, but we put it out to avoid having only one huge class

  private final HealpixNested hn;
  private final FillingCurve2D fc;

  HealpixNestedNeighbourSelector(HealpixNested healpixNested) {
    this.hn = healpixNested;
    this.fc = this.hn.fc;
  }

  private final void checkHashRange(long hash) {
    if (hash < 0 || hn.nHash <= hash) {
      throw new IllegalArgumentException("Hash value " + hash + " must be in [0, " + hn.nHash + "[");
    }
  }
  
  @Override
  public long neighbour(long hash, MainWind direction) {
    checkHashRange(hash);
    final HashParts hashParts = this.hn.decodeRegularHash(hash);
    return neighbour(hashParts, direction);
  }
  
  private long neighbour(final HashParts hashParts, MainWind direction) {
    return neighbour(hashParts.baseCellHash(), hashParts.iInBaseCell(), hashParts.jInBaseCell(), direction);
  }
  
  private long neighbour(int baseCellHash, int iInBaseCell, int jInBaseCell, MainWind direction) {
    final BaseHash baseHash = BaseHashes.get(baseCellHash);
    iInBaseCell += direction.getOffsetSE();
    jInBaseCell += direction.getOffsetSW();
    final MainWind neigBaseCellDirection = getNeighbourBaseCellDirection(iInBaseCell, jInBaseCell);
    baseCellHash = iInBaseCell;
    iInBaseCell = baseHash.pickRightIndexOnNeighbourSouthToEastAxis(neigBaseCellDirection,
        iInBaseCell, jInBaseCell, this.hn.nsideRemainderMask);
    jInBaseCell = baseHash.pickRightIndexOnNeighbourSouthToWestAxis(neigBaseCellDirection,
        baseCellHash, jInBaseCell, this.hn.nsideRemainderMask);
    baseCellHash = baseHash.getNeighbour(neigBaseCellDirection);
    return hn.encodeRegularHash(baseCellHash, iInBaseCell, jInBaseCell);
  }

  private MainWind getNeighbourBaseCellDirection(final int iInBaseCell, final int jInBaseCell) {
    final int neigBaseCellOffsetSE = neighbourBaseCellOffset(iInBaseCell);
    final int neigBaseCellOffsetSW = neighbourBaseCellOffset(jInBaseCell);
    return MainWind.getFromOffset(neigBaseCellOffsetSE, neigBaseCellOffsetSW);
  }

  private int neighbourBaseCellOffset(final int cooInBaseCell) {
    final int off = (cooInBaseCell >> 31) | (cooInBaseCell >> this.hn.depth);
    assert (cooInBaseCell == -1       && off == -1)
        || (cooInBaseCell == hn.nside && off ==  1)
        || (0 <= cooInBaseCell && cooInBaseCell < hn.nside && off == 0);
    return off;
  }
  
  @Override
  public NeighbourList neighbours(final long hash) {
    final NeighbourList result = new NeighbourList(hn.depth);
    this.neighbours(hash, result);
    return result;
  }
  
  @Override
  public void neighbours(long hash, final FlatHashList result) {
    checkHashRange(hash);
    result.clear();
    final HashBits hBits = this.hn.pullBitsApart(hash);
    if (isInBasePixelBorderFromBits(hBits.iInD0hBits, hBits.jInD0hBits)) {
      edgePixelNeighbours(hash, result);
    } else {
      innerPixelNeighbours(hBits.d0hBits, hBits.iInD0hBits, hBits.jInD0hBits, result);
    }
  }
  
  @Override
  public void neighbours(long hash, NeighbourList result) {
    checkHashRange(hash);
    result.clear();
    final HashBits hBits = this.hn.pullBitsApart(hash);
    if (isInBasePixelBorderFromBits(hBits.iInD0hBits, hBits.jInD0hBits)) {
      edgePixelNeighbours(hash, result);
    } else {
      innerPixelNeighbours(hBits.d0hBits, hBits.iInD0hBits, hBits.jInD0hBits, result);
    }
  }

  @Override
  public NeighbourList neighbours(final long hash, final EnumSet<MainWind> directions) {
    final NeighbourList result = new NeighbourList(hn.depth);
    this.neighbours(hash, directions, result);
    return result;
  }
  
  @Override
  public void neighbours(final long hash, final EnumSet<MainWind> directions, final NeighbourList result) {
    checkHashRange(hash);
    final HashBits hBits = this.hn.pullBitsApart(hash);
    result.clear();
    if (isInBasePixelBorderFromBits(hBits.iInD0hBits, hBits.jInD0hBits)) {
      edgePixelNeighbours(hash, directions, result);
    } else {
      innerPixelNeighbours(hBits.d0hBits, hBits.iInD0hBits, hBits.jInD0hBits, directions, result);
    }
  }

  /*private boolean isInBasePixelBorderFromCoo(final int iInBasePixel, final int jInBasePixel) {
    return 0 == iInBasePixel || iInBasePixel == this.hn.nsideRemainderMask
        || 0 == jInBasePixel || jInBasePixel == this.hn.nsideRemainderMask;
  }*/

  private boolean isInBasePixelBorderFromBits(
      final long iInBasePixelBits, final long jInBasePixelBits) {
    return 0 == iInBasePixelBits || iInBasePixelBits == this.hn.xMask
        || 0 == jInBasePixelBits || jInBasePixelBits == this.hn.yMask;
  }

  private void innerPixelNeighbours(long d0hBits, long iBits, long jBits,
      final FlatHashList result) {
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
    final long jm1Bits = ij & this.hn.yMask;
    final long im1Bits = ij & this.hn.xMask;
    // Compute i+1 and j+1 bits.
    // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
    // result in making fewer operations
    ij = this.fc.ij2hash(i + 1, j + 1);
    final long jp1Bits = ij & this.hn.yMask;
    final long ip1Bits = ij & this.hn.xMask;
    result.put(bits2hash(d0hBits, im1Bits, jm1Bits));
    result.put(bits2hash(d0hBits,   iBits, jm1Bits));
    result.put(bits2hash(d0hBits, ip1Bits, jm1Bits));
    result.put(bits2hash(d0hBits, im1Bits,   jBits));
    result.put(bits2hash(d0hBits, ip1Bits,   jBits));
    result.put(bits2hash(d0hBits, im1Bits, jp1Bits));
    result.put(bits2hash(d0hBits,   iBits, jp1Bits));
    result.put(bits2hash(d0hBits, ip1Bits, jp1Bits));
  }
  
  private void innerPixelNeighbours(long d0hBits, long iBits, long jBits,
      final NeighbourList result) {
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
    final long jm1Bits = ij & this.hn.yMask;
    final long im1Bits = ij & this.hn.xMask;
    // Compute i+1 and j+1 bits.
    // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
    // result in making fewer operations
    ij = this.fc.ij2hash(i + 1, j + 1);
    final long jp1Bits = ij & this.hn.yMask;
    final long ip1Bits = ij & this.hn.xMask;
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
    final long jm1Bits = ij & this.hn.yMask;
    final long im1Bits = ij & this.hn.xMask;
    // Compute i+1 and j+1 bits.
    // Again, depending to the FillingCurve implementation, calling 2 x fc.i02hash(...) twice could
    // result in making fewer operations
    ij = this.fc.ij2hash(i + 1, j + 1);
    final long jp1Bits = ij & this.hn.yMask;
    final long ip1Bits = ij & this.hn.xMask;
    for (final MainWind mw : directions) {
      final long nHash = hn.bits2hash(d0hBits,
          mw.pickRightSouthToEastLongValue(im1Bits, iBits, ip1Bits),
          mw.pickRightSouthToWestLongValue(jm1Bits, jBits, jp1Bits));
      addNeighbour(nHash, MainWind.SE, result);
    }
  }
  
  private void edgePixelNeighbours(long hash, final FlatHashList result) {
    final HashParts hashParts = this.hn.decodeRegularHash(hash);
    edgePixelNeighbours(hashParts, result);
  }
  
  private void edgePixelNeighbours(long hash, final NeighbourList result) {
    final HashParts hashParts = this.hn.decodeRegularHash(hash);
    edgePixelNeighbours(hashParts, result);
  }
  
  private void edgePixelNeighbours(final HashParts hashParts, final FlatHashList result) {
    // Could have simply been edgePixelNeighbours(hash, EnumSet.allOf(MainWind.class) result)
    // but we prefered to unroll the for loop.
    addNeighbourIfExists(neighbour(hashParts, MainWind.S ), result);
    addNeighbour        (neighbour(hashParts, MainWind.SE), result);
    addNeighbourIfExists(neighbour(hashParts, MainWind.E ), result);
    addNeighbour        (neighbour(hashParts, MainWind.SW), result);
    addNeighbour        (neighbour(hashParts, MainWind.NE), result);
    addNeighbourIfExists(neighbour(hashParts, MainWind.W ), result);
    addNeighbour        (neighbour(hashParts, MainWind.NW), result);
    addNeighbourIfExists(neighbour(hashParts, MainWind.N ), result);
  }
  
  private void edgePixelNeighbours(final HashParts hashParts, final NeighbourList result) {
    // Could have simply been edgePixelNeighbours(hash, EnumSet.allOf(MainWind.class) result)
    // but we prefered to unroll the for loop.
    addNeighbourIfExists(neighbour(hashParts, MainWind.S ), MainWind.S , result);
    addNeighbour        (neighbour(hashParts, MainWind.SE), MainWind.SE, result);
    addNeighbourIfExists(neighbour(hashParts, MainWind.E ), MainWind.E , result);
    addNeighbour        (neighbour(hashParts, MainWind.SW), MainWind.SW, result);
    addNeighbour        (neighbour(hashParts, MainWind.NE), MainWind.NE, result);
    addNeighbourIfExists(neighbour(hashParts, MainWind.W ), MainWind.W , result);
    addNeighbour        (neighbour(hashParts, MainWind.NW), MainWind.NW, result);
    addNeighbourIfExists(neighbour(hashParts, MainWind.N ), MainWind.N , result);
  }
  
  private void edgePixelNeighbours(long hash, final EnumSet<MainWind> directions,
      final NeighbourList neighbours) {
    final HashParts hashParts = this.hn.decodeRegularHash(hash);
    for (final MainWind direction : directions) {
      final long neigHash = neighbour(hashParts, direction);
      addNeighbourIfExists(neigHash, direction, neighbours);
    }
  }

  private static void addNeighbour(final long neighbourHash, final MainWind direction,
      final NeighbourList neighbours) {
    neighbours.put(neighbourHash, direction);
    /*neighbours.indirections[direction.getIndex()] = neighbours.size;
    neighbours.neighbours[neighbours.size++] = neighbourHash;*/
    
  }
  private static void addNeighbourIfExists(final long neighbourHash, final MainWind direction,
      final NeighbourList neighbours) { // Branch version
    if (neighbourHash >= 0) {
      addNeighbour(neighbourHash, direction, neighbours);
    }
  }
  private static void addNeighbour(final long neighbourHash, final FlatHashList result) {
    result.put(neighbourHash);
    
  }
  private static void addNeighbourIfExists(final long neighbourHash, final FlatHashList result) {
    if (neighbourHash >= 0) {
      addNeighbour(neighbourHash, result);
    }
  }
  /**
   * May be faster than {@link #addNeighbourIfExists} since according to 
   * http://ithare.com/infographics-operation-costs-in-cpu-clock-cycles/
   * a "Wrong" branch of "if" (branch miss-prediction) cost 10-20 CPU cycle.
   * Here, no "if" but 5 additional operations, each costing <1 cycles.
   * And also 5 + 4 useless operation when if is wong
   * @param neighbourHash
   * @param direction
   * @param neighbours
   */
  /*private void addNeighbourIfExistsV2(final long neighbourHash, final MainWind direction,
      final NeighbourList neighbours) { // Branch free version
    final int bits0or1 = (int) (neighbourHash >> 63);
    assert (neighbourHash == -1 && bits0or1 == -1)  || (neighbourHash >= 0 && bits0or1 == 0);
    
    neighbours.indirections[direction.getIndex()] = (neighbours.size | bits0or1);
    assert (neighbours.size | -1) == -1 && (neighbours.size | 0) == neighbours.size; 
    
    neighbours.neighbours[neighbours.size++] = neighbourHash;
    neighbours.size += bits0or1;
  }*/

  // DEALS WITH INTERNAL EDGES
  
  @Override
  public FlatHashList internalEdges(final long hash, final int deltaDepth) {
    assert 1 < deltaDepth && deltaDepth < 30;
    final int n = ((1 << deltaDepth) - 1) << 2; // 4 * nside - 4
    final FlatHashList result = new FlatHashList(this.hn.depth + deltaDepth, n);
    internalEdges(hash, deltaDepth, result);
    return result;
  }

  @Override
  public FlatHashList sortedInternalEdges(final long hash, final int deltaDepth) {
    assert 1 < deltaDepth && deltaDepth < 30;
    final int n = ((1 << deltaDepth) - 1) << 1; // 4 * nside - 4
    final FlatHashList result = new FlatHashList(this.hn.depth + deltaDepth, n);
    sortedInternalEdges(hash, deltaDepth, result);
    return result;
  }

  @Override
  public FlatHashList sortedInternalEdge(final long hash, final int deltaDepth,
      final Ordinal direction) {
    final FlatHashList result = new FlatHashList(this.hn.depth + deltaDepth, 1 << deltaDepth);
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
    final FlatHashList result = new FlatHashList(this.hn.depth + deltaDepth, 1 << deltaDepth);
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
    final FlatHashList result = new FlatHashList(this.hn.depth + deltaDepth, 1 << deltaDepth);
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
    final FlatHashList result = new FlatHashList(this.hn.depth + deltaDepth, 1 << deltaDepth);
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
    final FlatHashList result = new FlatHashList(this.hn.depth + deltaDepth, 1 << deltaDepth);
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

  /**
   * Write in the given {@link FlatHashList} the hashes corresponding to the
   * internal bounds of the hash at the hash depth + deltaDepth.
   * The first quarter contains the southeast border (the z-order curve x-axis
   * with y = 0), the second quarter contains the northeast border
   * (the z-order y-axis with x = xmax - 1), the third quarter contains the
   * northwest border (the z-order curve x-axis with y = ymax - 1) and the
   * forth quarter contains the southwest border (the y-axis with x = 0).
   * The hashes are ordered consecutively, starting from the south (x=0, y=0)
   * cell in the anti-clokwise direction.
   * @param hash the hash for which we look for the internal bounds
   * @param toEdgeDeltaDepth difference between the depth of the edge cells
   * and the depth of the given cell.
   * @param result the list used to store the result (its content, if any,
   *  is overwritten).
   */
  @Override
  public void internalEdges(long hash, final int toEdgeDeltaDepth,
      final FlatHashList result) {
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

  /**
   * This is for instance useful if you want to access the data sequentially
   * from a HDD.
   *
   * @param hash
   * @param deltaDepth
   * @param result
   */
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

  // DEALS WITH EXTERNAL EDGES

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
  
  void appendSortedInternalEdgeElement(final long hash, final int toEdgeDeltaDepth,
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
  
  @Override
  public FlatHashList externalEdges(long hash, int deltaDepth) {
    final FlatHashList result =  new FlatHashList(this.hn.depth + deltaDepth, 4 + (4 << deltaDepth));
    externalEdges(hash, deltaDepth, result);
    return result;
  }
  
  @Override
  public void externalEdges(long hash, int toEdgeDeltaDepth, final FlatHashList result) {
    externalEdges(hash, toEdgeDeltaDepth, result, false);
  }
  
  @Override
  public FlatHashList sortedExternalEdges(long hash, int deltaDepth) {
    final FlatHashList result =  new FlatHashList(this.hn.depth + deltaDepth, 4 + (4 << deltaDepth));
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
    final NeighbourList neihbours = new NeighbourList(hn.depth);
    final HashBits hBits = this.hn.pullBitsApart(hash);
    if (isInBasePixelBorderFromBits(hBits.iInD0hBits, hBits.jInD0hBits)) {
      // Not easy: opposite directions depends on base cell neighbours
      final HashParts hashParts = this.hn.decodeRegularHash(hash);
      final BaseHash baseHash = BaseHashes.get(hashParts.baseCellHash());
      edgePixelNeighbours(hashParts, neihbours);
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
      innerPixelNeighbours(hBits.d0hBits, hBits.iInD0hBits, hBits.jInD0hBits, neihbours);
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

  // UTILIY METHODS
  
  private static long xMask(final int depth) {
    return 0x5555555555555555L >>> (64 - (depth << 1));
  }

  private static long yMask(final int depth) {
    return 0xAAAAAAAAAAAAAAAAL >>> (64 - (depth << 1));
  }
  
  private static long xyMask(final int depth) {
    return (-1L) >>> (64 - (depth << 1));
  }

}
