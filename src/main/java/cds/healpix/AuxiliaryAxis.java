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
 * Define an auxiliary axis, i.e. a 3rd dimension to be taken into account when computing
 * an HEALPix index.<br>
 * To compute a hash using the auxiliary axis, the values must be projected along the z-axis
 * ranging from 0 (inclusive) to 1 (exclusive).<br>
 * The maximum depth when using an auxiliary axis is 19 (4 + 19 * 3 = 61 bits).
 * So the maximum number of intervals is 2^19 = 524288.<br>
 * At a given {@code depth}, the [0.0, 1.0[ interval is divided into {@code NSIDE} (=2^depth) segments.<br>
 * The resolution at a given {@code depth} is thus abs(unproject(1) - unproject(0)) / {@code NSIDE}.<br>
 * If we want to represent temporal data, we can represent:
 * <ul>
 *   <li>about 1 years with intervals being time ranges of 1 min;</li>
 *   <li>60.68 years with intervals being time ranges of 1 hour.</li>
 * </ul>
 * We recall that when a double is divided by a power of two, there is no loss in precision since
 * the mantissa remains unchanged (only the exponent part changes).<br>
 * 
 * Example: 
 *   depth 12 =&gt; NSIDE = 4096 =&gt; angular  resolution = 51.5 arcsec
 *            =&gt;  4096 / 365 = 11.222 years with a resolution  = 1 day
 *
 * @author F.-X. Pineau
 *
 */
public interface AuxiliaryAxis {

  /**
   * Returns the name of the auxiliary axis.
   * @return the name of the auxiliary axis.
   */
  String name();
  
  /**
   * Project (or normalize) the given auxiliary axis value into the [0.0, 1.0[ range.
   * @param x the auxiliary axis value to be projected (normalized).
   * @return the projected (normalized) value on the auxiliary axis, in [0.0, 1.0[.
   */
  double project(double x);
  
  /**
   * Unproject (denormalize) the value from the [0.0, 1.0[ range.
   * @param y the projected (normalized) value on the auxiliary axis, in [0.0, 1.0[
   * @return the unprojected (denormalize) value on the auxiliary axis.
   */
  double unproject(double y);
  
}
