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

import cds.healpix.AbstractHealpixNestedMOC.CurrentCellAccessor;

/**
 * WARNING: IN DEVELOPPEMTENT, NOT YET VISIBLE IN THE JAVADOC
 * TODO: CONTINUE DEV
 * 
 * @author F.-X. Pineau
 *
 */
final class HealpixNestedMOC extends AbstractHealpixNestedMOC<CurrentCellAccessor> {

  private HealpixNestedMOC(final int mocDepth, final long[] mocCells) {
    this(mocDepth, mocCells, mocCells.length);
  }

  private HealpixNestedMOC(final int mocDepth, final long[] mocCells, final int toIndex) {
    super(mocDepth, mocCells, toIndex);
  }

  public HealpixNestedMOC toDeeperDepth(final int newDepth) {
    if (newDepth <= this.depthMax) {
      throw new IllegalArgumentException("The given depth must be higher than the depth max of the MOC");
    }
    final int twiceDepthDiff = (newDepth - this.depthMax) << 1; 
    final long[] mocNew = new long[this.to];
    for (int i = 0; i < this.to; i++) {
      long c = this.cells[i];
      mocNew[i] = c << twiceDepthDiff;
    }
    return new HealpixNestedMOC(newDepth, mocNew);
  }
  
  public HealpixNestedMOC toLowerDepth(final int newDepth) {
    if (newDepth >= this.depthMax) {
      throw new IllegalArgumentException("The given depth lust be lower than the depth max of the MOC");
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
          (currHashAtNewDepth = getHashFromDepth(raw, depth) >> (depth - newDepth))) {
        mocNew[iNew++] = (prevHashAtNewDepth << 1) | 1L;
        prevHashAtNewDepth = currHashAtNewDepth;
      }
    }
    return new HealpixNestedMOC(newDepth, Arrays.copyOf(mocNew, iNew));
  }

  /**
   * The depth of the returned MOC is the larger depth between this MOC and the given MOC.
   * @param rightMoc
   * @return
   */
  /*public HealpixNestedMOC innerJoin(final HealpixNestedMOC rightMoc) { // = intersect
    final int depthMax = Math.max(this.depthMax, rightMoc.depthMax);
    final long[] inter = new long[Math.max(this.size(), rightMoc.size())];
    int iInter = 0;
    final Iterator<CurrentCellAccessor> leftIt = this.iterator();
    final Iterator<CurrentCellAccessor> rightIt = rightMoc.iterator();
    CurrentCellAccessor currLeft, currRight;
    if (!leftIt.hasNext() || !rightIt.hasNext()) {
      return new HealpixNestedMOC(depthMax, new long[]{});
    } else {
      currLeft = leftIt.next();
      currRight = rightIt.next();
    }
    for (;;) {
      if (currLeft.getDepth() < currRight.getDepth()) {
        //if ()
      } else if (currLeft.getDepth() > currRight.getDepth()) {
        break;
      } else {
        assert currLeft.getDepth() == currRight.getDepth();
        
      }
    }
    
    return null; //new HealpixNestedMOC(newDepth, Arrays.copyOf(mocNew, iNew));
  }*/
  
  
  
  /**
   * Create a MOC considering that the given array is already a MOC: i.e. it is sorted (ASC order)
   * and do not contains duplicate or small cells included into larger one's.
   * WARNING: the array is used internally, so it must not be modified by an external reference!
   * use {@code Arrays.copy()} is you are not sure! 
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedMOC createUnsafe(final int mocDepth, final long[] mocCells) {
    return new HealpixNestedMOC(mocDepth, mocCells);
  }

  /**
   * Same as {@link #createUnsafe(int, long[])} except that not we do not use the full array.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @param toIndex the index of the last element (exclusive) to be considered in the moc
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedMOC createUnsafe(final int mocDepth, final long[] mocCells, final int toIndex) {
    return new HealpixNestedMOC(mocDepth, mocCells);
  }

  /**
   * Same a {@link #createUnsafe(int, long[])} except that the properties (array sorted,
   * no duplicates, no cell included into an other one) is checked.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedMOC createCheck(final int mocDepth, final long[] mocCells) {
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
  public static HealpixNestedMOC createCheck(final int mocDepth, final long[] mocCells, final int toIndex) {
    checkIsSortedNoDuplicateNoOverlap(mocCells);
    return createUnsafe(mocDepth, mocCells, toIndex);
  }

  public static HealpixNestedBMOC create(int mocDepth, long[] mocCells) {
    // also make a version taking a sortedIterator on hash of depth < deptMax
    throw new Error("Method not yet implemented!");
    // sort 
    // remove duplicate
    // if ambiguity (flag), generate an error
  }

  @Override
  public Iterator<CurrentCellAccessor> iterator() {
    return new Iter() {
      @Override
      protected CurrentCellAccessor returnThis() {
        return this;
      }
      @Override
      protected void calledInNext(int i) { }
      @Override
      public void remove() {
       throw new UnsupportedOperationException();
      }
    };
  }

}
