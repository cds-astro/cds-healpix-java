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
import java.util.NoSuchElementException;

import cds.healpix.CompassPoint.MainWind;

/**
 * Class storing the hash values of the neighbour cells of a cell of given hash, together with the
 * direction of each neighbour with respect to the central cell.
 * 
 * @author F.-X. Pineau
 *
 */
public final class NeighbourList implements ListOfHash {
  
  // We do not use a EnumMap because of the primitive 'long':
  // we want to avoid autoboxing!

  private final int depth;
  private final long[] neighbours = new long[MainWind.size()];
  private int size;
  private final MainWind[] neighMainWinds = new MainWind[MainWind.size()];
  private final int[] indirections = new int[MainWind.size()]; // contains for each MainWind the corresponding index in the neighbours array
  
  
  public NeighbourList(final int depth) {
    this.depth = depth;
    clear();
  }

  public void clear() {
    Arrays.fill(this.indirections, -1);
    Arrays.fill(this.neighMainWinds, null);
    this.size = 0;
  }
  
  public boolean contains(final MainWind mainWind) {
    final int i = this.indirections[mainWind.getIndex()];
    return i != -1;
  }

  public long get(final MainWind mainWind) {
    final int i = this.indirections[mainWind.getIndex()];
    return i == -1 ? i : this.neighbours[i];
  }  
  
  void put(final long neighbourHash, final MainWind mainWind) {
    this.neighbours[size] = neighbourHash;
    this.neighMainWinds[size] = mainWind;
    this.indirections[mainWind.getIndex()] = size;
    ++size;
  }
  
  public void sortByHashAsc() {
    if (size > 1) {
      final int sizem1 = size - 1;
      // Sort hashes and main winds
      for (int from = 0; from < sizem1; from++) {
        int iMin = indexOfMinValue(neighbours, from, this.size);
        if(iMin != from) {
          // Swap hashes at indices 'from' and 'iMin'
          long tmpH = neighbours[from];
          neighbours[from] = neighbours[iMin];
          neighbours[iMin] = tmpH;
          // Swap mainwinds at indices 'from' and 'iMin'
          MainWind tmpMW = neighMainWinds[from];
          neighMainWinds[from] = neighMainWinds[iMin];
          neighMainWinds[iMin] = tmpMW;
        }
      }
      // Update indirections
      for (int i = 0; i < size; i++) {
        final MainWind mw = neighMainWinds[i];
        indirections[mw.getIndex()] = i;
      }
    }
  }
  
  private static final int indexOfMinValue(long[] array, int from, int to) {
    int iMin = from;
    long vMin = array[from];
    for (int i = from + 1; i < to; i++) {
      if (array[i] < vMin) {
        vMin = array[i];
        iMin = i;
      }
    }
    return iMin;
  }

  @Override
  public int size() {
    return this.size;
  }

  @Override
  public long get(final int i) {
    checkIndex(i);
    return this.neighbours[i];
  }
  
  public MainWind getDirection(final int i) {
    checkIndex(i);
    return neighMainWinds[i];
  }

  private final void checkIndex(final int i) {
    if (i >= this.size) {
      throw new IndexOutOfBoundsException("i: " + i + " ; size: " + size);
    }
  } 
  
  @Override
  public void arraycopy(int srcPos, long[] dest, int destPos, int length) {
    if (srcPos + length > this.size) {
      throw new IndexOutOfBoundsException("srcPos + length > source size");
    }
    System.arraycopy(this.neighbours, srcPos, dest, destPos, length);
  }

  @Override
  public FlatHashIterator iterator() {
    return new FlatHashIterator() {
      private int i = 0;
      @Override
      public int depth() {
        return depth;
      }
      @Override
      public boolean hasNext() {
        return this.i < size;
      }
      @Override
      public long next() {
        if (this.i >= size) {
          throw new NoSuchElementException();
        }
        return neighbours[this.i++];
      }
    };
  }
}
