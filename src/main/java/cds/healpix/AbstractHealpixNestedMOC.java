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

import java.util.Iterator;

abstract class AbstractHealpixNestedMOC<E extends AbstractHealpixNestedMOC.CurrentCellAccessor> implements Iterable<E>  {
  //Internally, a long is formed of bits:
  // BBBBxx...xxS00...00F if depth < depthMax
  // BBBBxx...xxxx...xxSF if depth = dephtMax
  // With
  //  -  B: the 4 bits coding the base hash [0- 11]
  //  - xx: the 2 bits of level x
  //  -  S: the sentinel bit coding the depth
  //  - 00: if (depth != depthMax) those bits are unused bits
  //  -  F; the flag bit (0: partial, 1: full)
  // Properties:
  // - a smaller depth cell is not necessarily after or before a  deeper depth cell it contains !!
  //   => to build the MOC, sort first by ignoring the Sentinel  bit and the flag (and,
  //      in case of equality, sort in DESC on the sentinel bit only):  easier to then remove cells
  //      included into a larger one)
  // WARNING: at depth 29, all 64 bits of a long are used. WE MUST CHECK WHAT IT MEANS WHEN SORTING BY NATURAL ORDER!

  protected final int depthMax;
  protected final long[] cells;
  protected final int to;

  protected AbstractHealpixNestedMOC(final int mocDepth, final long[] mocCells) {
    this(mocDepth, mocCells, mocCells.length);
  }

  protected AbstractHealpixNestedMOC(final int mocDepth, final long[] mocCells, final int toIndex) {
    checkDepth(mocDepth);
    this.depthMax = mocDepth;
    this.cells = mocCells;
    this.to = toIndex;
  }
  
  protected static final void checkDepth(final int depth) {
    if (depth >= Healpix.DEPTH_MAX) {
      throw new IllegalArgumentException("Depth larger than or equal to " +  Healpix.DEPTH_MAX
          + " not uspported.");
    }
  }
  
  protected static long rmSentinel(long mocEncodedHash) {
    return turnOffRightmost1Bit(mocEncodedHash);
  }

  protected static long getHash(final long mocEncodedHash, final int depthMax) {
    return getHashFromDepth(mocEncodedHash, getDepth(mocEncodedHash, depthMax));
  }

  protected static long getHashFromDepth(final long mocEncodedHash, final int depth) {
    return mocEncodedHash >>> (1 + (depth << 1));
  }

  protected static int getDepth(long mocEncodedHash, final int depthMax) {
    return depthMax - (Long.numberOfTrailingZeros(mocEncodedHash) >> 1); // number of leading zero / 2
  }

  protected static final boolean aIsALargerCellThanB(long mocEncodedHashA, long mocEncodedHashB) {
    return aHasSmallerDepthThanB(mocEncodedHashA, mocEncodedHashB);
  }

  protected static final boolean aHasSmallerDepthThanB(long mocEncodedHashA, long mocEncodedHashB) {
    return isolateRightmost1Bit(mocEncodedHashA) < isolateRightmost1Bit(mocEncodedHashB);
  }

  protected static final long turnOffRightmost1Bit(final long x) {
    return x & (x - 1);
  }

  protected static final long isolateRightmost1Bit(final long x) {
    return x & -x;
  }
  
  public static final long encodeHash4MOC(int depth, long hash, int depthMax) {
    // Set the sentinel bit
    hash <<= 1;
    hash |= 1L;
    // Shit according to the depth and add space for the flag bit
    hash <<= (1 + ((depthMax - depth) << 1));
    return hash;
  }
  
  protected static final boolean contains(final long largerCellMocEncodedHashNoSentinel,
      final long smallerCellMocEncodedHashNoSentinek) {
    return (largerCellMocEncodedHashNoSentinel & smallerCellMocEncodedHashNoSentinek)
        == largerCellMocEncodedHashNoSentinel;
  }

  protected static void checkIsSortedNoDuplicateNoOverlap(final long[] mocCells) {
    if (mocCells.length > 0) {
      long prev = mocCells[0];
      long prevNoSentinel = rmSentinel(prev);
      for (int i = 1; i < mocCells.length; i++) {
        final long curr = mocCells[i];
        final long currNoSentinel = rmSentinel(curr);
        // Test is sorted
        if (prev >= curr) {
          throw new IllegalArgumentException("Not valid MOC: not sorted");
        }
        // Test no overlap (if ordered, smaller depth always before the elements it contains)
        long aAndB = prevNoSentinel & currNoSentinel;
        if (aAndB == prevNoSentinel || aAndB == currNoSentinel) {
          throw new IllegalArgumentException("Not valid MOC: elem " + curr + " overlap with " + prev);
        }
        // loop
        prev = curr;
        prevNoSentinel = currNoSentinel;
      }
    }
  }
  
  /**
   * Returns the maximal depth of this MOC.
   * @return the maximal depth of this MOC.
   */
  public final int getDepthMax() {
    return this.depthMax;
  }
  
  /**
   * Returns the number of elements the moc contains, i.e. the number of cells at various depth.
   * @return the number of elements the moc contains, i.e. the number of cells at various depth.
   */
  public final int size() {
    return this.to;
  }
  
  /**
   * Returns the number of cells at depth {@link #getDepthMax()} the moc contains, i.e.
   * the sum for each cell of the number of cells at depth {@link #getDepthMax()}.
   * @return the number of cells at depth {@link #getDepthMax()} the moc contains, i.e.
   * the sum for each cell of the number of cells at depth {@link #getDepthMax()}.
   */
  public final long computeDeepSize() {
    long tot = 0;
    for (int i = 0; i < this.to; i++) {
      final int depth = getDepth(this.cells[i], this.depthMax);
      tot += Healpix.nsideSquare(depth);
    }
    return tot;
  }
  
  public static interface CurrentCellAccessor {
    long getMOCEncodedHash();
    int getDepth();
    long getHash();
  }

  protected abstract class Iter implements Iterator<E>, CurrentCellAccessor {

    private int i = 0;

    private long h = -1;
    private int d = 0;

    protected abstract E returnThis();
    protected abstract void calledInNext(int i);
    
    @Override
    public final long getMOCEncodedHash() {
      return this.h;
    }
    @Override
    public final int getDepth() {
      return depthMax - this.d;
    }
    @Override
    public final long getHash() {
      return this.h >> (1 + (this.d << 1)); // remove 1 bits per depth + 1 sentinel bit
    }
    @Override
    public final boolean hasNext() {
      return this.i < to;
    }
    @Override
    public final E next() {
      this.h = cells[this.i];
      this.d = Long.numberOfTrailingZeros(h) >> 1;
      this.calledInNext(i);
      this.i++;
      return this.returnThis();
    }
    @Override
    public String toString() {
      return String.format("d: %d; h: %d; raw: %d", getDepth(), getHash(), getMOCEncodedHash());
    }
  }

}
