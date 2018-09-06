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

import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.abs;
import static cds.healpix.common.math.Math.sqrt;

/**
 * Class defines (and storing the vertices of) a polygon on the unit sphere.
 * 
 * @author F.-X. Pineau
 *
 */
public final class Polygon {

  /**
   * Defines the method used to know if the south pole is in the polygon or in its complement.
   *
   * @author F.-X. Pineau
   *
   */
  public enum ContainsSouthPoleComputer {
    /**
     * We explicitly tell that the south pole is inside the polygon.
     */
    PROVIDED_TRUE() {
      @Override
      boolean containsSouthPole(final CooXYZ[] vertices, final Vect3D[] normals) {
        return true;
      }
    },
    /**
     * We explicilty tell that the south pole is NOT inside the polygon.
     */
    PROVIDED_FALSE() {
      @Override
      boolean containsSouthPole(final CooXYZ[] vertices, final Vect3D[] normals) {
        return true;
      }
    },
    /**
     * Consider the south pole to be inside the polygon if the sum of consecutive longitude
     * differences equals 2pi and if there is more vertices in the south hemisphere than in the 
     * north hemisphere.
     */
    BASIC() {
      @Override
      boolean containsSouthPole(final CooXYZ[] vertices, final Vect3D[] normals) {
        int nVertexInSouthHemiSphere = 0;
        double sumDeltaLon = 0;
        for (int i = 0, j = vertices.length - 1; i < vertices.length; j = i++) {
          final double deltaLon = vertices[i].lon() - vertices[j].lon();
          final double absDeltaLon = abs(deltaLon);
          if (absDeltaLon <= PI) {
            sumDeltaLon += deltaLon;
          } else if (deltaLon > 0) {
            sumDeltaLon -= 2 * PI - absDeltaLon;
          } else {
            assert deltaLon < 0;
            sumDeltaLon += 2 * PI - absDeltaLon;
          }
          if (vertices[i].lat() < 0) {
            nVertexInSouthHemiSphere++;
          }
        }
        return abs(sumDeltaLon) > PI // sumDeltaLon = 0 or -2PI or 2PI
            && (nVertexInSouthHemiSphere << 1) > vertices.length; // more vertices in south that in north
      }
    },
    /**
     * Consider that the gravity center of the first 3 non-aligned vertices is inside the polygon
     * if it is on the left of the fisrt two edges, then test is the south pole is inside or outside
     * the polygon.
     */
    STD_FXP() {
      @Override
      boolean containsSouthPole(final CooXYZ[] vertices, final Vect3D[] normals) {
        final Vect3D firstEdgeNormal = CooXYZ.crossProd(vertices[0], vertices[1]);
        int i = 1;
        while (abs(vertices[++i].scalarProd(firstEdgeNormal)) < 1e-10); // bug for a circular polygon made of very small segments
        CooXYZ triangleCenter = CooXYZ.normalizedSum(vertices[0], vertices[1], vertices[i]);
        return triangleCenter.scalarProd(firstEdgeNormal) > 0
            ^ oddNumberOfIntersectionGoingSouth(
                new CooXYZ(triangleCenter.x(), triangleCenter.y(), triangleCenter.z()),
                vertices, normals);
      }
    },
    /**
     * Consider that the inside of the polygon is always located on the left part of each edge
     * (thus, the inside become the outside when considering the vertices list in reverse order),
     * then test is the south pole is inside or outside the polygon.
     * WARNING: Valid only for non self-intersecting polygons.
     * WARNING: NOT IMPLENTED YET (I HAVE TO SEE HOW TO IMPLEMENT IT FOR CONCAVE POLYGONS).
     */
    STD_IVOA() {
      @Override
      boolean containsSouthPole(final CooXYZ[] vertices, final Vect3D[] normals) {
        throw new Error("Not implemented yet!");
      }
    };
    abstract boolean containsSouthPole(final CooXYZ[] vertices, final Vect3D[] normals);
  }
  
  
  /** Coordinates of the North pole. */
  private static final CooXYZ NORTH_POLE = new CooXYZ(PI, 0.5 * PI);

  /** Array of the polygon vertices. */
  private final CooXYZ[] vertices;

  /** Cross products of consecutive vertices. */
  private final Vect3D[] vectProds;

  /** Tells if the polygon contains one of the two poles. */
  private final boolean containsSouthPole;

  /**
   * Create a new polygon from the given list of vertices, using the
   * {@link ContainsSouthPoleComputer#BASIC} method to define its inside and outside.
   * @param polyVertices vertices defining the polygon
   */
  public Polygon(final CooXYZ[] polyVertices) {
    this(polyVertices, ContainsSouthPoleComputer.BASIC);
  }
  
  /**
   * Create a new polygon from the given list of vertices, using the
   * given method to define its inside and outside.
   * @param polyVertices vertices defining the polygon
   * @param cspc method used to defined the inseide and the outside of the polygon
   */
  public Polygon(final CooXYZ[] polyVertices, ContainsSouthPoleComputer cspc) {
    if (polyVertices.length < 3) {
      throw new IllegalArgumentException("A polygon must have"
          + " a minimum of 3 vertices!");
    }
    this.vertices = polyVertices;
    this.vectProds = new Vect3D[this.vertices.length];
    for (int i = 0, j = this.vertices.length - 1; i < this.vertices.length; j = i++) {
      final Vect3D xprod = CooXYZ.crossProd(this.vertices[i], this.vertices[j]);
      if (NORTH_POLE.scalarProd(xprod) < 0) {
        this.vectProds[i] = xprod.opposite();
      } else {
        this.vectProds[i] = xprod;
      }
    }
    this.containsSouthPole = cspc.containsSouthPole(polyVertices, this.vectProds);
  }
  
  /**
   * Returns the number of vertices the polygon contains.
   * @return the number of vertices the polygon contains.
   */
  public int nVertices() {
    return this.vertices.length;
  }

  /**
   * Returns the vertex located at the given index in the polygin vertex list.
   * @param vertexIndex index of the vertex we want to access (in [0, nVertices[
   * @return the vertex located at the given index in the polygin vertex list.
   */
  public CooXYZ vertex(final int vertexIndex) {
    return this.vertices[vertexIndex];
  }

  /**
   * Returns {@code true} if the polygon contain the point p.
   * @param p point to test
   * @return {@code true} if the polygon contain the point p.
   */
  public boolean contains(final CooXYZ p) {
    return this.containsSouthPole ^ oddNumberOfIntersectionGoingSouth(p);
  }

  /**
   * Returns {@code true} if there is a odd number of intersection going from the given point
   * to the south pole (at the given piont constant longitude).
   * @param p point to test
   * @return {@code true} if there is a odd number of intersection going from the given point
   * to the south pole (at the given piont constant longitude).
   */
  private boolean oddNumberOfIntersectionGoingSouth(final CooXYZ p) {
    boolean c = false;
    for (int i = 0, j = this.vertices.length - 1; i < this.vertices.length; j = i++) {
      if (isInLonRange(p,  this.vertices[j],  this.vertices[i])
          && crossPlaneGoingSouth(p, this.vectProds[i])) {
        c = !c;
      }
      /*// Commented the time to verify results are correct
       if ((   (abs(this.vertices[i].lon() - this.vertices[j].lon()) > PI)
            ^ (this.vertices[i].lon() > p.lon() != this.vertices[j].lon() > p.lon())
          ) && p.scalarProd(this.vectProds[i]) > 0) {
        c = !c;
      }*/
    }
    return c;
  }
  
  private static boolean oddNumberOfIntersectionGoingSouth(
      final CooXYZ p, final CooXYZ[] vertices, final Vect3D[] normals) {
    boolean c = false;
    for (int i = 0, j = vertices.length - 1; i < vertices.length; j = i++) {
      if (isInLonRange(p,  vertices[j],  vertices[i])
          && crossPlaneGoingSouth(p, normals[i])) {
        c = !c;
      }
    }
    return c;
  }
  
  /**
   * Returns {@code true} if the given point p longitude is between the given vertices v1 and v2
   * longitude range
   */
  private static boolean isInLonRange(final CooXYZ p, final CooXYZ v1, final CooXYZ v2) {
    // Case do not cross lon == 0
    //    0 < v2.lon - v1.lon <= PI && p > v1.lon && p < v2.lon
    // || 0 < v1.lon - v2.lon <= PI && p > v2.lon && p < v1.lon
    // Cross lon == 0
    // || v2.lon - v1.lon > PI && (p < v1.lon || p > v2.lon)
    // || v1.lon - v2.lon > PI && (p < v2.lon || p > v1.lon)
    return (abs(v2.lon() - v1.lon()) > PI) ^ (v2.lon() > p.lon() != v1.lon() > p.lon());
  }
  /**
   * Returns {@code true} if the line at constant (x, y) and decreasing z going from the given point
   * toward south intersect the plane of given normal vector. The normal vector must have a positive
   * z coordinate (=> mist be in the north hemisphere)
   */
  private static boolean crossPlaneGoingSouth(final CooXYZ p,
      final Vect3D planeNormalDirInNorthHemisphere) {
    return p.scalarProd(planeNormalDirInNorthHemisphere) > 0;
  }
  
  
  /**
   * Returns {@code true} if an edge of the polygone intersects the line defined by the two given
   * points.
   * @param a first segment point
   * @param b second segment point
   * @return {@code true} if an edge of the polygone intersects the line defined by the two given
   * points.
   */
  public boolean intersectSegAB(CooXYZ a, CooXYZ b) {
    double ua, ub;
    // Ensure a < b in longitude
    if (a.lon() > b.lon()) {
      CooXYZ swp = a;
      a = b;
      b = swp;
    }
    for (int i = 0, j = this.vertices.length - 1; i < this.vertices.length; j = i++) {
      CooXYZ pA = this.vertices[i];
      CooXYZ pB = this.vertices[j];
      // Ensures pA < pB in longitude
      if (pA.lon() > pB.lon()) {
        CooXYZ swp = pA;
        pA = pB;
        pB = swp;
      }
      if (segmentsAreOverlappingInLon(a, b, pA, pB) 
          && polygonEdgeIntersectsGreatCircle(ua = a.scalarProd(this.vectProds[i]), ub = b.scalarProd(this.vectProds[i]))
          && intersectPointInPolygonSegment(a, b, pA, pB, this.vectProds[i], ua, ub)) {
        return true;
      }
    }
    return false;
  }
  
  private static boolean segmentsAreOverlappingInLon(
      final CooXYZ a, final CooXYZ b, final CooXYZ pA, final CooXYZ pB) {
    return (pB.lon() - pA.lon() > PI)
        ^ (a.lon() <= pB.lon() && b.lon() >= pA.lon() && (b.lon() - a.lon()) <= PI);
  }
  
  /**
   * Tells if segment from vector a to vector b intersect the plane of the great circle defined
   * by its normal vector N.
   * @param aDotProdEdgeN the dot product of vector a with the great circle normal vector N.
   * @param bDotProdEdgeN the dot product of vector b with the great circle normal vector N.
   * @return {@code true} if vector a and b are in opposite part of the plane having for normal
   * vector vector N.
   */
  private static boolean polygonEdgeIntersectsGreatCircle(
      final double aDotProdEdgeN, final double bDotProdEdgeN) {
    return (aDotProdEdgeN > 0) != (bDotProdEdgeN > 0);
  }
  
  /**
   * Tells if the intersetion line (i) between the two plane defined by vector a, b and pA, pB
   * respectively in inside the zone [pA, pB].
   * @param a
   * @param b
   * @param pA
   * @param pB
   * @param aDotProdEdgeN
   * @param bDotProdEdgeN
   * @return
   */
  private static boolean intersectPointInPolygonSegment(
      final CooXYZ a, final CooXYZ b, final CooXYZ pA, final CooXYZ pB, final Vect3D paCrossProdPb,
      final double aDotProdEdgeN, final double bDotProdEdgeN) {
    final Vect3D intersect = normalizedIntersectPoint(a, b, paCrossProdPb, aDotProdEdgeN, bDotProdEdgeN);
    final double papb = pA.scalarProd(pB);
    return abs(pA.scalarProd(intersect)) > papb && abs(pB.scalarProd(intersect)) > papb;
  }
  
  private static Vect3D normalizedIntersectPoint(
      final CooXYZ a, final CooXYZ b, final Vect3D paCrossProdPb,
      final double aDotProdEdgeN, final double bDotProdEdgeN) {
    // We note u = a x b
    // Intersection vector i defined by
    // i = (pA x pB) x (a x b)
    //   = u x (a x b)
    //   = (u.b)a - (u.a)b
    // i = i / ||i|| 
    final double x = bDotProdEdgeN * a.x() - aDotProdEdgeN * b.x();
    final double y = bDotProdEdgeN * a.y() - aDotProdEdgeN * b.y();
    final double z = bDotProdEdgeN * a.z() - aDotProdEdgeN * b.z();
    final double norm = sqrt(x * x + y * y + z * z);
    return new Vect3D(x / norm, y / norm, z / norm);
  }
  
}
