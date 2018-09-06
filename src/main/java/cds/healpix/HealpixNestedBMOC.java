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
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * A BMOC is a MOC storing for each element a binary value telling if the cell is PARTIALY or
 * FULLY covered by a region.<br>
 * Internally, a long is made of bits:<br>
 * <b> BBBBxx...xxS00...00F if depth &lt; depthMax </b><br>
 * <b> BBBBxx...xxxx...xxSF if depth = dephtMax </b><br>
 * With:<br>
 * <ul>
 *   <li> B: the 4 bits coding the base hash [0- 11]</li>
 *   <li>xx: the 2 bits of level x</li>
 *   <li> S: the sentinel bit coding the depth</li>
 *   <li>00: if (depth != depthMax) those bits are unused bits</li>
 *   <li> F; the flag bit (0: partial, 1: full)</li>
 * </ul>
 * WARNING: not tested yet at depth 29, may not work because of Java signed long
 *          (all 64 bits of a long are used: 62 + 1 sentinel + 1 flag)
 * 
 * @author F.-X. Pineau
 *
 */
public final class HealpixNestedBMOC implements Iterable<HealpixNestedBMOC.CurrentValueAccessor> {

  // a smaller depth cell is not necessarily after or before a  deeper depth cell it contains !!
  // => to build the MOC, sort first by ignoring the Sentinel bit and the flag (and,
  //    in case of equality, sort in DESC on the sentinel bit only):  easier to then remove cells
  //    included into a larger one)
  // WARNING: at depth 29, all 64 bits of a long are used. WE MUST CHECK WHAT IT MEANS WHEN SORTING BY NATURAL ORDER!

  private static final long MASK_EXCEPT_LSB = -1L << 1; // <=> ~1L

  /**
   * Tells if an element of the MOC is covered or not by a region and, if covered, it tells is
   * the element is partially or fully covered. 
   * @author F.-X. Pineau
   *
   */
  public static enum Status {
    /** The cell does not cover the region. */
    EMPTY,
    /** The cell is partially covered by the region. */
    PARTIAL,
    /** The cell is fully covered by the region. */
    FULL
  }

  /**
   * Provides informations on the current element of the BMOC while iterating over it.
   * @author F.-X. Pineau
   *
   */
  public static interface CurrentValueAccessor {
    /**
     * 
     * @return the current BMOC value, encoded in a specific way
     */
    long getRawValue();
    /**
     * Returns the depth of the BMOC element this cursor is pointing at.
     * @return the depth of the BMOC element this cursor is pointing at.
     */
    int getDepth();
    /**
     * Returns then hash value of the BMOC element this cursor is pointing at.
     * @return then hash value of the BMOC element this cursor is pointing at.
     */
    long getHash();
    /**
     * Returns {@code true} if the status the BMOC element this cursor is pointing at is
     * {@link Status#FULL}.
     * @return {@code true} if the status the BMOC {@code true} if the status the BMOC element this cursor is pointing at is
     * {@link Status#FULL} element this cursor is pointing at is
     * {@link Status#FULL}.
     */
    boolean isFull();
  }

  private final class Iter implements Iterator<CurrentValueAccessor>, CurrentValueAccessor {

    private int i = 0;

    private long h = -1;
    private int d = 0;

    @Override
    public long getRawValue() {
      return h;
    }
    @Override
    public int getDepth() {
      return depthMax - d;
    }
    @Override
    public long getHash() {
      return h >> (2 + (d << 1)); // remove 2 bits per depth + 1 sentinel bit + 1 flag bit 
    }
    @Override
    public boolean isFull() {
      return HealpixNestedBMOC.isFull(h);
    }
    @Override
    public boolean hasNext() {
      return this.i < to;
    }
    @Override
    public CurrentValueAccessor next() {
      this.h = cells[this.i++];
      this.d = Long.numberOfTrailingZeros(h >> 1) >> 1; // First remove the flag bit, then /= 2
      return this;
    }
    @Override
    public String toString() {
      return String.format("d: %d; h: %d; full: %s; raw: %d", getDepth(), getHash(), isFull(), getRawValue());
    }
    @Override
    public void remove() {
     throw new UnsupportedOperationException();
    }
  }

  private final int depthMax;
  private final long[] cells;
  private final int to;

  private HealpixNestedBMOC(final int mocDepth, final long[] mocCells) {
    this(mocDepth, mocCells, mocCells.length);
  }
  
  private HealpixNestedBMOC(final int mocDepth, final long[] mocCells, final int toIndex) {
    checkDepth(mocDepth);
    this.depthMax = mocDepth;
    this.cells = mocCells;
    this.to = toIndex;
  }
  
  /**
   * Returns a new BMOC having a deeper depth.
   * @param newDepth the depeht of the wanted new BMOC
   * @return a new BMOC having a deeper depth.
   */
  public HealpixNestedBMOC toDeeperDepth(final int newDepth) {
    if (newDepth <= this.depthMax) {
      throw new IllegalArgumentException("The given depth must be higher than the depth max of the MOC");
    }
    final int twiceDepthDiff = (newDepth - this.depthMax) << 1; 
    final long[] mocNew = new long[this.to];
    for (int i = 0; i < this.to; i++) {
      long c = this.cells[i];
      mocNew[i] = ((c & MASK_EXCEPT_LSB) << twiceDepthDiff) | (c & 1L);
    }
    return new HealpixNestedBMOC(newDepth, mocNew);
  }
  
  public HealpixNestedBMOC toLowerDepth(final int newDepth) { // toDeeperDepth: only change depthMax!!
    if (newDepth >= this.depthMax) {
      throw new IllegalArgumentException("The given depth must be lower than the depth max of the MOC");
    }
    final long[] mocNew = new long[this.to];
    int iNew = 0;
    long prevHashAtNewDepth = -1, currHashAtNewDepth;
    for (int i = 0; i < this.to; i++) {
      final long raw = this.cells[i];
      final int depth = getDepth(raw, this.depthMax);
      if (depth <= newDepth) {
        mocNew[iNew++] = raw;
      } else if (prevHashAtNewDepth != 
          (currHashAtNewDepth = getHashFromDepthDiff(raw, this.depthMax - depth) >> (depth - newDepth))) {
          mocNew[iNew++] = (prevHashAtNewDepth << 2) | 2L; // sentinel bit + flag = 0
          prevHashAtNewDepth = currHashAtNewDepth;
      }
    }
    return new HealpixNestedBMOC(newDepth, Arrays.copyOf(mocNew, iNew));
  }
 
  private void checkDepth(final int depth) {
    if (depth >= Healpix.DEPTH_MAX) {
      throw new IllegalArgumentException("Depth larger than or equal to " +  Healpix.DEPTH_MAX
          + " not uspported.");
    }
  }
  
  
  /**
   * Returns the number of elements the moc contains, i.e. the number of cells at various depth.
   * @return the number of elements the moc contains, i.e. the number of cells at various depth.
   */
  public int size() {
    return this.to;
  }
  
  /**
   * Returns the number of cells at depth {@link #getDepthMax()} the moc contains, i.e.
   * the sum for each cell of the number of cells at depth {@link #getDepthMax()}.
   * @return the number of cells at depth {@link #getDepthMax()} the moc contains, i.e.
   * the sum for each cell of the number of cells at depth {@link #getDepthMax()}.
   */
  public long computeDeepSize() {
    long tot = 0;
    for (int i = 0; i < this.to; i++) {
      final int depth = getDepth(this.cells[i], this.depthMax);
      tot += Healpix.nsideSquare(this.depthMax - depth);
    }
    return tot;
  }
  
  /**
   * Create a MOC considering that the given array is already a MOC: i.e. it is sorted (ASC order)
   * and do not contains duplicate or small cells included into larger one's.
   * WARNING: the array is used internally, so it must not be modified by an external reference!
   * use {@code Arrays.copy()} is you are not sure! 
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedBMOC createUnsafe(final int mocDepth, final long[] mocCells) {
    return new HealpixNestedBMOC(mocDepth, mocCells);
  }

  /**
   * Same as {@link #createUnsafe(int, long[])} except that not we do not use the full array.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @param toIndex the index of the last element (exclusive) to be considered in the moc
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedBMOC createUnsafe(final int mocDepth, final long[] mocCells, final int toIndex) {
    return new HealpixNestedBMOC(mocDepth, mocCells, toIndex);
  }
  
  /**
   * Same a {@link #createUnsafe(int, long[])} except that the properties (array sorted,
   * no duplicates, no cell included into an other one) is checked.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedBMOC createCheck(final int mocDepth, final long[] mocCells) {
    checkIsSortedNoDuplicateNoOverlap(mocCells);
    return createUnsafe(mocDepth, mocCells);
  }

  /**
   * Same as {@link #createCheck(int, long[])} except that not we do not use the full array.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @param toIndex the index of the last element (exclusive) to be considered in the moc
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedBMOC createCheck(final int mocDepth, final long[] mocCells, final int toIndex) {
    checkIsSortedNoDuplicateNoOverlap(mocCells);
    return createUnsafe(mocDepth, mocCells, toIndex);
  }
  
  
  public static HealpixNestedBMOC createPacking(int mocDepth, long[] mocCells) {
    return createPacking(mocDepth, mocCells, mocCells.length);
  }
  
  /**
   * We assume here that the given array is ordered, that no element overlaps another one, but that
   * the MOC is not normalized, i.e. a large cell may be splitted in 4 sub-cell (each sub-cell
   * possibly splitted in its 4 sub-cells recursively). 
   * @param mocDepth depth of the moc
   * @param mocCells ordered list of cells
   * @param toIndex index of the last cell to be read in the given array of cells
   * @return a new moc from the input parameters, packing if necessary.
   */
  public static HealpixNestedBMOC createPacking(int mocDepth, long[] mocCells, final int toIndex) {
    // On-place pack
    int prevToIndex = 0;
    int currToIndex = toIndex;
    while (prevToIndex != currToIndex) { // changes occurs
      prevToIndex = currToIndex;
      int iPrevMoc = 0, iCurrMoc = 0;
      while (iPrevMoc < prevToIndex) {
        long currCell = mocCells[iPrevMoc++];
        int currCellDepth = getDepth(currCell, mocDepth);
        long currCellHash = getHashFromDepthDiff(currCell, mocDepth - currCellDepth);
        // Look for the first cell of the larger cell (depth - 1)  (=> 2 last bits = 00), the cell must be FULL
        while (iPrevMoc < prevToIndex &&
            (currCellDepth == 0 || isPartial(currCell) || isNotFirstCellOfLargerCell(currCellHash))) {
          if (iCurrMoc != iPrevMoc) {
            mocCells[iCurrMoc++] = currCell;
          }
          currCell = mocCells[iPrevMoc++];
          currCellDepth = getDepth(currCell, mocDepth);
          currCellHash = getHashFromDepthDiff(currCell, mocDepth - currCellDepth);
        }
        // Look at the 3 sibling
        if (iPrevMoc + 2 < prevToIndex
            && mocCells[iPrevMoc + 0] == buildValue(currCellDepth, currCellHash | 1, true, mocDepth)
            && mocCells[iPrevMoc + 1] == buildValue(currCellDepth, currCellHash | 2, true, mocDepth)
            && mocCells[iPrevMoc + 2] == buildValue(currCellDepth, currCellHash | 3, true, mocDepth)) {
          mocCells[iCurrMoc++] = buildValue(currCellDepth - 1, currCellHash >> 2, true, mocDepth);
          iPrevMoc += 3;
        } else if (iCurrMoc != iPrevMoc) {
          mocCells[iCurrMoc++] = currCell;
        }
      }
      currToIndex = iCurrMoc;
    }
    // We may find a better algorithm doing a sngle pass on the input MOC
    // Here the number of passes max = mocDepth - smallestDepthOfACellInOutputMoc
    return createUnsafe(mocDepth, mocCells, currToIndex);
  }
  
  private static final boolean isNotFirstCellOfLargerCell(final long hash) {
    return (hash & 3L) != 0L;
  }
  
  /*public static HealpixNestedBMOC create(int mocDepth, long[] mocCells) {
    // also make a version taking a sortedIterator on hash of depth < deptMax
    
    throw new Error("Method not yet implemented!");
    // sort 
    // remove duplicate
    // if ambiguity (flag), generate an error
  }*/

  private static void checkIsSortedNoDuplicateNoOverlap(final long[] mocCells) {
    if (mocCells.length > 0) {
      long prev = rmFlag(mocCells[0]);
      long prevNoSentinel = rmSentinelNoFlag(prev);
      for (int i = 1; i < mocCells.length; i++) {
        final long curr = rmFlag(mocCells[i]);
        final long currNoSentinel = rmSentinelNoFlag(curr);
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

  private static long rmFlag(long rawVal) {
    return rawVal >> 1;
  }

  private static long setFlagTo0(long rawVal) {
    return rawVal & MASK_EXCEPT_LSB;
  }

  private static long rmSentinelNoFlag(long valNoFlag) { // or flag removed, or set to 0
    assert (valNoFlag & 1L) == 0L;
    return turnOffRightmost1Bit(valNoFlag);
  }

  private static long getHash(final long rawVal, final int depthMax) {
    return getHashFromDepthDiff(rawVal, depthMax - getDepth(rawVal, depthMax));
  }
  
  private static long getHashFromDepthDiff(final long rawVal, final int depthDiff) {
      return rawVal >> (2 + (depthDiff << 1));
  }
  
  private static int getDepth(long rawVal, final int depthMax) {
    return getDepthNoFlag(rawVal >> 1, depthMax); // Remove the flag bit, then /= 2
  }
  
  private static int getDepthNoFlag(long valNoFlag, final int depthMax) { // or flag removed, or set to 0
    // ssert (valNoFlag & 1L) == 0L;
    return depthMax - (Long.numberOfTrailingZeros(valNoFlag) >> 1); // number of leading zero / 2 (works both with no flag or with flag set to 0)
  }

  private static boolean isPartial(long mocEncodedHash) {
    return (mocEncodedHash & 1L) == 0L; 
  }
  
  private static boolean isFull(long mocEncodedHash) {
    return (mocEncodedHash & 1L) == 1L; 
  }
  
  /**
   * Returns the BMOC deeper depth.
   * @return the BMOC deeper depth.
   */
  public int getDepthMax() {
    return this.depthMax;
  }

  /**
   * Returns the status of the given hash at the given depth
   * @param depth depth of the hash we want the status
   * @param hash  hash for which we want the status
   * @return the status of the given hash at the given depth
   */
  public Status statut(int depth, long hash) {
    checkDepth(depth);
    // assert cells is sorted!!!
    long h = buildValue(depth, hash, false, this.depthMax);
    long hNoSentinel = rmSentinelNoFlag(h);
    if (h < this.cells[0]) {
      final long h0 =  this.cells[0];
      final long h0NoFlag = setFlagTo0(h0);
      if (h == h0NoFlag) {
        return Status.FULL; // Was smaller only because the flag bit is set on cells[0]
      }
      final long h0NoSentinelNoFlag = rmSentinelNoFlag(h0NoFlag);
      if (aIsALargerCellThanB(h, h0NoFlag)) {
        return contains(hNoSentinel, h0NoSentinelNoFlag) ? Status.PARTIAL : Status.EMPTY;
      } else {
        return contains(h0NoSentinelNoFlag, hNoSentinel) ? selfStatus(h0, h0NoFlag) : Status.EMPTY; 
      }
    } else if (h >  this.cells[this.to - 1]) {
      final long hn =  this.cells[this.to - 1];
      final long hnNoFlag = setFlagTo0(hn);
      final long hnNoSentinelNoFlag = rmSentinelNoFlag(hnNoFlag);
      if (aIsALargerCellThanB(h, hnNoFlag)) {
        return contains(hNoSentinel, hnNoSentinelNoFlag) ? Status.PARTIAL : Status.EMPTY;
      } else {
        return contains(hnNoSentinelNoFlag, hNoSentinel) ? selfStatus(hn, hnNoFlag) : Status.EMPTY; 
      }
    } else {
      int index = Arrays.binarySearch(this.cells, 0, this.to, h);
      if(index >= 0) {
        return Status.PARTIAL; // since we found the same value, with flag set to 0
      }
      // look at the larger elem)
      long hp1 = this.cells[-index - 1];
      long hp1NoFlag = setFlagTo0(hp1);
      if (h == hp1NoFlag) {
        return Status.FULL;
      }
      final long hp1NoSentinelNoFlag = rmSentinelNoFlag(hp1NoFlag);
      if (aIsALargerCellThanB(h, hp1NoFlag)) {
        if (contains(hNoSentinel, hp1NoSentinelNoFlag)) {
          return Status.PARTIAL;
        }
      } else if (contains(hp1NoSentinelNoFlag, hNoSentinel) ) {
        return selfStatus(hp1, hp1NoFlag); 
      }
      // look at the smaller elem
      long hm1 = this.cells[-index - 2];
      long hm1NoFlag = setFlagTo0(hm1);
      if (h == hm1NoFlag) {
        return Status.FULL;
      }
      final long hm1NoSentinelNoFlag = rmSentinelNoFlag(hm1NoFlag);
      if (aIsALargerCellThanB(h, hm1NoFlag)) {
        if (contains(hNoSentinel, hm1NoSentinelNoFlag)) {
          return Status.PARTIAL;
        }
      } else if (contains(hm1NoSentinelNoFlag, hNoSentinel) ) {
        return selfStatus(hm1, hm1NoFlag); 
      }
      return Status.EMPTY;
    }
  }

  private static Status selfStatus(long val) {
    return (val & 1L) == 0L ? Status.PARTIAL : Status.FULL;
  }
  
  private static Status selfStatus(long val, long valFlagSetTo0) {
    return val == valFlagSetTo0 ? Status.PARTIAL : Status.FULL;
  }
  
  /**
   * Creates a specific hash encoding the hash value at the given depth, the depth and the status
   * flag.
   * @param depth depth of the hash.
   * @param hash value of the hash.
   * @param isFull flag telling if the cell is fully covered or not.
   * @param depthMax maximum depth of the MOC the value will belong to.
   * @return a BMOC encoded hash value with its depth and status flag.
   */
  public static final long buildValue(int depth, long hash, boolean isFull, int depthMax) {
    // Set the sentinel bit
    hash <<= 1;
    hash |= 1L;
    // Shit according to the depth and add space for the flag bit
    hash <<= (1 + ((depthMax - depth) << 1));
    // Set the flag bit
    hash |= (isFull ? 1L : 0L);
    return hash;
  }

  private static final boolean aIsALargerCellThanB(long aNoFlag, long bNoFlag) {
    return aHasSmallerDepthThanB(aNoFlag, bNoFlag);
  }
  
  private static final boolean aHasSmallerDepthThanB(long aNoFlag, long bNoFlag) {
    return isolateRightmost1Bit(aNoFlag) < isolateRightmost1Bit(bNoFlag);
  }

  private static final long turnOffRightmost1Bit(final long x) {
    return x & (x - 1);
  }

  private static final long isolateRightmost1Bit(final long x) {
    return x & -x;
  }
  
  private static final boolean contains(final long largerCellHashNoFlagNoSentinel,
      final long smallerCellHashNoFlagNoSentinek) {
   return (largerCellHashNoFlagNoSentinel & smallerCellHashNoFlagNoSentinek)
        == largerCellHashNoFlagNoSentinel;
  }

  @Override
  public Iterator<CurrentValueAccessor> iterator() {
    return new Iter();
  }
  
  /**
   * Returns an iterator on all the cells in the BMOC at the maximum depth.
   * @return an iterator on all the cells in the BMOC at the maximum depth.
   */
  public FlatHashIterator flatHashIterator() {
    return new FlatHashIterator() {
      private int i = 0;
      private long currVal = to == 0 ? -1 : init();
      private long currMax;
      
      private long init() {
        final long mocEncodedHash = cells[i++];
        final int depth = getDepth(mocEncodedHash, depthMax);
        final int depthDiff = depthMax - depth;
        final long hash = getHashFromDepthDiff(mocEncodedHash, depthDiff);
        final int twiceDepthDiff = depthDiff << 1;
        this.currVal = hash << twiceDepthDiff;
        this.currMax = currVal | ((1L << twiceDepthDiff) - 1L);
        return currVal;
      }
      
      private void internalNext() {
        if (this.currVal < this.currMax) {
          ++this.currVal;
        } else if (this.i < to) {
          init();
        } else {
          this.currVal = -1;
        }
      }
      
      @Override
      public int depth() {
        return depthMax;
      }
      @Override
      public boolean hasNext() {
        return this.currVal != -1;
      }
      @Override
      public long next() {
        final long res = this.currVal;
        if (res == -1) {
          throw new NoSuchElementException();
        }
        internalNext();
        return res;
      }
    };
  }

}
