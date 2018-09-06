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

import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.sqrt;
import static cds.healpix.common.math.Math.atan2;

/**
 * Three dimensional vector.
 * @author F.-X. Pineau
 *
 */
final class Vect3D {
    
    /** X coordinate. */
    private final double x;
    
    /** Y coordinate. */
    private final double y;
    
    /** Z coordinate. */
    private final double z;
    
    /**
     * Constructor for the Cartesian coordiantes.
     * @param x first Cartesian coordinate
     * @param y second Cartesian coordinate
     * @param z thrid Cartesian coordinate
     */
    public Vect3D(final double x, final double y, final double z) {
        this.x = x;
        this.y = y; 
        this.z = z;
    }
     
    /**
     * Constructor from the spherical coordinates on the unit sphere.
     * @param lonRad
     * @param latRad
     */
    public Vect3D(final double lonRad, final double latRad) {
        final double cosDec = Math.cos(latRad);
        this.x = cosDec * Math.cos(lonRad);
        this.y = cosDec * Math.sin(lonRad);
        this.z = Math.sin(latRad);
    }
    
    /**
     * Getter
     * @return the x cartesian coordinate
     */
    public final double x() { return this.x; }
    /**
     * Getter
     * @return the y cartesian coordinate
     */
    public final double y() { return this.y; }
    /**
     * Getter
     * @return the z cartesian coordinate
     */
    public final double z() { return this.z; }
    
    public Vect3D clone(){
        return new Vect3D(this.x, this.y, this.z);
    }
    
    /**
     * Returns the longitude coordinate from the vector, in radians.
     * @return the longitude coordinate from the vector, in radians.
     */
    public double toLon(){
        final double tmp = Math.atan2(this.y, this.x);
        return tmp < 0 ? 2 * PI + tmp : tmp ;
    }
    
    /**
     * Returns the latitude coordinates from the vector, in radians.
     * @return the latitude coordinates from the vector, in radians.
     */
    public double toLat(){
        return atan2(this.z, sqrt(this.x * this.x + this.y * this.y));
    }
    
    /**
     * Returns a vector having the opposite coordinates of this vector.
     * @return a vector having the opposite coordinates of this vector.
     */
    public Vect3D opposite() {
        return new Vect3D(-this.x, -this.y, -this.z);
    }
    
    /**
     * Returns the norm of this vector.
     * @return the norm of this vector.
     */
    public double norm(){
        return sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }
    
    /**
     * Returns the vector having the same direction than this vector,
     * but normalized (norm = 1).
     * @return the vector having the same direction than this vector,
     * but normalized (norm = 1).
     */
    public final Vect3D normalized(){
        final double norm = this.norm();
        return new Vect3D(this.x / norm, this.y / norm, this.z / norm);
    }
    
    /**
     * Returns the scalar product of v1 by v2.
     * @param v1 vector 1
     * @param v2 vector 2
     * @return  the scalar product of v1 by v2.
     */
    public static double scalarProd(final Vect3D v1, final Vect3D v2) {
        return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
    }
    
    /**
     * Returns the cross product of v1 by v2.
     * @param v1 vector 1
     * @param v2 vector 2
     * @return the cross product of v1 by v2.
     */
    public static Vect3D crossProd(final Vect3D v1, final Vect3D v2) {
        return new Vect3D((v1.y * v2.z) - (v1.z * v2.y)
                ,(v1.z * v2.x) - (v1.x * v2.z)
                ,(v1.x * v2.y) - (v1.y * v2.x));
    }
    
    @Override
    public final String toString() {
        return String.format(Locale.US, "(%+011.9f, %+011.9f, %+011.9f)",
                this.x, this.y, this.z); 
    }
}
