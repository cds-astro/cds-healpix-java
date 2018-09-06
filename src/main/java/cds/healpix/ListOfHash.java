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

public interface ListOfHash extends FlatHashIterable {

  /**
   * Returns the current number of elements in the flat list.
   * @return the current number of elements in the flat list.
   */
  int size();
 
  /**
   * Returns the element at the given index {@code i} in the list.
   * @param i index in the list of the element to be returned.
   * @return the element at the given index {@code i} in the list.
   */
  long get(final int i);
  
  /**
   * Similar to {@link System#arraycopy(Object, int, Object, int, int)} except that the source
   * is the internal list.
   * @param srcPos see {@link System#arraycopy(Object, int, Object, int, int)}
   * @param dest see {@link System#arraycopy(Object, int, Object, int, int)}
   * @param destPos see {@link System#arraycopy(Object, int, Object, int, int)}
   * @param length see {@link System#arraycopy(Object, int, Object, int, int)}
   */
  void arraycopy(int srcPos, long[] dest, int destPos, int length);
  
}
