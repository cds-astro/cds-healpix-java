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

/**
 * Define an iterator over the primitive type 'long'.
 * We do this because Iterators in Java can have Object generic types, not primitive types, and
 * we want to avoid auto-boxing for performances reasons.
 * 
 * Additionally, we use this interface to iterate over hashes of a same depth given by 
 * the {@link #depth()} method.
 *
 * @author F.-X. Pineau
 *
 */
public interface FlatHashIterator extends HierarchyItem {

  /**
   * Returns {@code true} if the iteration has more elements. (In other words, returns {@code true}
   * if {@link #next()} would return an element rather than throwing an exception.)
   * @return {@code true} if the iteration has more elements
   */
  boolean hasNext();
  
  /**
   * Returns the next element in the iteration.
   * @return the next element in the iteration.
   */
  long next();
  
}
