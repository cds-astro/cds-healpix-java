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

/**
 * This class is used when the number of Hash returned y a function is not necessarily
 * known in advance but has an upper bound.
 * It is basically a read-only array list on long primitives.
 * 
 * @author F.-X. Pineau
 *
 */
public final class FlatHashList implements ListOfHash {

  // Used internally as a DataStructure and externally as an Object.
  
  final int depth;
  /** A simple array of long. */
  final long[] hList;
  /** Number of elements in the list (can be smaller than the list length). */
  int size = 0;

  /**
   * The capacity of the flat list, i.e. its maximum size, i.e the length of
   * the internal array used to store the list elements.
   * @param depth depth of the hash in the list.
   * @param capacity maximum size the list.
   */
  public FlatHashList(final int depth, final int capacity) {
    this.depth = depth;
    this.hList = new long[capacity];
  }

  /**
   * Returns the maximum number of elements that can be stored in this list.
   * @return the maximum number of elements that can be stored in this list.
   */
  public int capacity() {
    return this.hList.length;
  }

  public void clear() {
    this.size = 0;
  }
  public FlatHashList put(long hash) {
    this.hList[this.size++] = hash;
    return this;
  }
  public FlatHashList put(ListOfHash hashes) {
    hashes.arraycopy(0, this.hList, this.size, hashes.size());
    this.size += hashes.size();
    return this;
  }


  public void sortByHashAsc() {
    Arrays.sort(this.hList, 0, this.size);
  }

  @Override
  public int size() {
    return this.size;
  }

  @Override
  public long get(final int i) {
    if (i >= this.size) { // if i < 0, IndexOutOfBoundsException also thrown
      throw new IndexOutOfBoundsException("i: " + i + " ; size: " + size);
    }
    return this.hList[i];
  }

  @Override
  public void arraycopy(int srcPos, long[] dest, int destPos, int length) {
    if (srcPos + length > this.size) {
      throw new IndexOutOfBoundsException("srcPos + length > source size");
    }
    System.arraycopy(this.hList, srcPos, dest, destPos, length);
  }

  @Override
  public FlatHashIterator iterator() {
    return new FlatHashIterator(){
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
        return hList[this.i++];
      }
    };
  }
  
  @Override
  public String toString() {
    return Arrays.toString(Arrays.copyOf(this.hList, this.size));
  }
}
