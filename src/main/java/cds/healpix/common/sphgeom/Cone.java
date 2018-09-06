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

package cds.healpix.common.sphgeom;

import java.util.Locale;

/**
 * Defines a cone on the unit sphere.
 * 
 * @author F.-X. Pineau
 *
 */
public final class Cone extends CooXYZ {

  private final double r;

  /**
   * Creator.
   * @param center cetner of the cone
   * @param radiusRad radius of the cone, in radians
   */
  public Cone(final CooXYZ center, final double radiusRad) {
    super(center.x(), center.y(), center.z());
    this.r = radiusRad;
  }

  /**
   * Creator.
   * @param lonRad longitude of the center of the cone, in radians
   * @param latRad latitude of the center of the cone, in radians
   * @param radiusRad  radius of the cone, in radians
   */
  public Cone(final double lonRad, final double latRad, final double radiusRad) {
    super(lonRad, latRad);
    this.r = radiusRad;
  }

  /**
   * Creator.
   * @param x first Cartesian coordinate of the center of the cone
   * @param y second Cartesian coordinate of the center of the cone
   * @param z third Cartesian coordinate of the center of the cone
   * @param radiusRad  radius of the cone, in radians
   */
  public Cone(final double x, final double y, final double z,
      final double radiusRad) {
    super(x, y, z);
    this.r = radiusRad;
  }

  /**
   * Returns the angle of the cone (the distance between its center and the its edge), in radians.
   * @return the angle of the cone (the distance between its center and the its edge), in radians.
   */
  public double radiusRad() {
    return this.r;
  }

  /**
   * Returns {@code true} if the given point is inside the cone.
   * @param coo position to be tested
   * @return {@code true} if the given point is inside the cone.
   */
  public boolean contains(final CooXYZ coo) {
    return spheDist(this, coo) <= this.r;
  }

  @Override
  public String toString() {
    return String.format(Locale.US,
        "Center (%.6f°, %.6f°); Radius: %.6f arcmin",
        Math.toDegrees(this.lon()), Math.toDegrees(this.lat()),
        Math.toDegrees(this.r) * 60);
  }
} 
