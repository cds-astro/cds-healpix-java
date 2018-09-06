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
 * The interface defines an item which is part of a Hierarchy.
 * By convention, the top layer (or root) of the hierarchy is at depth 0 and the lower
 * the layer is, the higher the depth is.
 *
 * @author F.-X. Pineau
 *
 */
public interface HierarchyItem {
    /**
     * Returns the depth of the item in a Hierarchy.
     * By convention, the top layer of the hierarchy is at depth 0 and the
     * lower the layer is, the higher the depth is.
     * @return the depth of the item in a Hierarchy.
     */
    int depth();
}
