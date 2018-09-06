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
import java.util.EnumMap;

import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.common.sphgeom.Cone;
import cds.healpix.common.sphgeom.CooXYZ;
import cds.healpix.common.sphgeom.Polygon;

import static cds.healpix.CompassPoint.Cardinal.E;
import static cds.healpix.CompassPoint.Cardinal.N;
import static cds.healpix.CompassPoint.Cardinal.S;
import static cds.healpix.CompassPoint.Cardinal.W;
import static cds.healpix.HealpixNestedBMOC.buildValue;
import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;
import static cds.healpix.Healpix.getBestStartingDepth;

public final class HealpixNestedPolygonComputer {

  private final int depthMax;
  private final HealpixNested hnDepthMax;
  private final HashComputer hcDepthMax;
  private final HashComputer[] hcs;
  private final VerticesAndPathComputer[] vpcs;
  private final NeighbourSelector[] neiSelect;

  private final EnumMap<Cardinal, double[]> vertices = new EnumMap<Cardinal, double[]>(Cardinal.class);
  private CooXYZ vertexN;
  private CooXYZ vertexE;
  private CooXYZ vertexS;
  private CooXYZ vertexW;
  private final FlatHashList neigs;
  
  private static interface AdditionalCheck {
    boolean isOk(Polygon poly, long hash);
  }
  private static final AdditionalCheck ALWAYS_OK = new AdditionalCheck() {
    @Override
    public boolean isOk(final Polygon poly, final long hash) {
      return true;
    }
  };
  private final AdditionalCheck centerInPoly;

  HealpixNestedPolygonComputer(final HealpixNested healpixNestedAtDepthMax) {
    this.depthMax = healpixNestedAtDepthMax.depth;
    this.hnDepthMax = healpixNestedAtDepthMax;
    this.hcDepthMax = this.hnDepthMax.newHashComputer();
    this.hcs = new HashComputer[this.depthMax + 1];
    this.vpcs = new VerticesAndPathComputer[this.depthMax + 1];
    this.neiSelect = new NeighbourSelector[this.depthMax + 1];
    this.neigs = new FlatHashList(this.depthMax, 12);
    this.vertices.put(N, new double[2]);
    this.vertices.put(E, new double[2]);
    this.vertices.put(S, new double[2]);
    this.vertices.put(W, new double[2]);
    this.centerInPoly = new AdditionalCheck() {
      private final VerticesAndPathComputer vpc = getVPC(depthMax);
      @Override
      public boolean isOk(final Polygon poly, final long hash) {
        final double[] ccoos = this.vpc.center(hash);
        final CooXYZ cellCenter = new CooXYZ(ccoos[LON_INDEX], ccoos[LAT_INDEX]);
        return poly.contains(cellCenter);
      }
    };
  }

  private HashComputer getHC(final int depth) {
    final int deltaDepth = this.depthMax - depth;
    HashComputer hc = this.hcs[deltaDepth];
    if (hc == null) {
      hc = Healpix.getNested(depth).newHashComputer();
      this.hcs[deltaDepth] = hc;
    }
    return hc;
  }

  private VerticesAndPathComputer getVPC(final int depth) {
    final int deltaDepth = this.depthMax - depth;
    VerticesAndPathComputer vpc = this.vpcs[deltaDepth];
    if (vpc == null) {
      vpc = Healpix.getNested(depth).newVerticesAndPathComputer();
      this.vpcs[deltaDepth] = vpc;
    }
    return vpc;
  }

  private NeighbourSelector getNeigSelect(final int depth) {
    final int deltaDepth = this.depthMax - depth;
    NeighbourSelector nei = this.neiSelect[deltaDepth];
    if (nei == null) {
      nei = Healpix.getNested(depth).newNeighbourSelector();
      this.neiSelect[deltaDepth] = nei;
    }
    return nei;
  }

  public HealpixNestedBMOC overlappingCells(double[][] vertices) {
    return overlapping(vertices, ALWAYS_OK);
  }

  public HealpixNestedBMOC overlappingCenters(double[][] vertices) {
   return overlapping(vertices, this.centerInPoly);
  }
  
  private HealpixNestedBMOC overlapping(double[][] vertices, final AdditionalCheck check) {
    final CooXYZ[] polyVertices = buildArrayOfCooXYZ(vertices);
    final Polygon poly = new Polygon(polyVertices);
    final Cone mec = CooXYZ.mec(polyVertices);
    int startingDepth = getBestStartingDepth(mec.radiusRad());
    // Compute smaller depth cells hashs
    if (startingDepth == -1) {
      neigs.clear();
      neigs.put(0).put(1).put(2).put(3).put(4).put(5).put(6).put(7).put(8).put(9).put(10).put(11);
      startingDepth = 0;
    } else {
      final long hashCenter = getHC(startingDepth).hash(mec.lon(), mec.lat());
      getNeigSelect(startingDepth).neighbours(hashCenter, neigs);
      neigs.put(hashCenter);
      neigs.sortByHashAsc();
    }
    // Compute and sort the list of cell containing at least one polygon vertex
    final long[] polyVerticesHash = buildArrayOfPolyVerticesHash(polyVertices);
    Arrays.sort(polyVerticesHash);
    // Build the list (removing duplicated) for all deltaDepth?
    final ArraysListOfLong moc = new ArraysListOfLong(10000);
    for (int i = 0 ; i < neigs.size; i++) {
      final long icell = neigs.get(i);
      buildMocRecursively(moc, startingDepth, icell, poly, polyVerticesHash, check);
    } 
    return HealpixNestedBMOC.createUnsafe(this.depthMax, moc.a, moc.size());
  }
  
  private void buildMocRecursively(final ArraysListOfLong moc, int depth, long hash,
      final Polygon poly, final long[] polyVerticesHash, final AdditionalCheck check) {
    int nVerticesInPoly;
    if (isInList(depth, hash, polyVerticesHash)
        || ((nVerticesInPoly = nVerticesInPoly(depth, hash, poly)) > 0 && nVerticesInPoly < 4)
        || hasIntersection(depth, hash, poly)) {
      if (depth == this.depthMax) {
        if (check.isOk(poly, hash)) {
          moc.add(buildValue(depth, hash, false, this.depthMax));
        }
      } else {
        hash <<= 2;
        depth++;
        buildMocRecursively(moc, depth,   hash, poly, polyVerticesHash, check);
        buildMocRecursively(moc, depth, ++hash, poly, polyVerticesHash, check);
        buildMocRecursively(moc, depth, ++hash, poly, polyVerticesHash, check);
        buildMocRecursively(moc, depth, ++hash, poly, polyVerticesHash, check);
      }
    } else if (nVerticesInPoly == 4) {
      moc.add(buildValue(depth, hash, true, this.depthMax));
    } // assert nVerticesInPoly == 0 && hasNoIntersection 
  }

  private boolean hasIntersection(int depth, long hash, final Polygon poly) {
   return  poly.intersectSegAB(this.vertexN, this.vertexE)
        || poly.intersectSegAB(this.vertexS, this.vertexE)
        || poly.intersectSegAB(this.vertexW, this.vertexN)
        || poly.intersectSegAB(this.vertexW, this.vertexS);
  }
  
  private int nVerticesInPoly(int depth, long hash, final Polygon poly) {
    final VerticesAndPathComputer vpc = this.getVPC(depth);
    vpc.vertices(hash, this.vertices);
    double[] vcoos = this.vertices.get(N);
    this.vertexN = new CooXYZ(vcoos[LON_INDEX], vcoos[LAT_INDEX]); // Costly: compute sin/cos
    vcoos = this.vertices.get(E);
    this.vertexE = new CooXYZ(vcoos[LON_INDEX], vcoos[LAT_INDEX]); // Costly: compute sin/cos
    vcoos = this.vertices.get(S);
    this.vertexS = new CooXYZ(vcoos[LON_INDEX], vcoos[LAT_INDEX]); // Costly: compute sin/cos
    vcoos = this.vertices.get(W);
    this.vertexW = new CooXYZ(vcoos[LON_INDEX], vcoos[LAT_INDEX]); // Costly: compute sin/cos
    int nIn = 0;
    if (poly.contains(this.vertexN)) { 
      nIn++; 
    }
    if (poly.contains(this.vertexE)) { 
      nIn++; 
    }
    if (poly.contains(this.vertexS)) { 
      nIn++; 
    }
    if (poly.contains(this.vertexW)) { 
      nIn++; 
    }
    return nIn;
  }

  private boolean isInList(int depth, long hash, final long[] polyVerticesHash) {
    final int twiceDeltaDepth = (this.depthMax - depth) << 1;
    long hashAtDepthMax = hash << twiceDeltaDepth;
    int i = Arrays.binarySearch(polyVerticesHash, hashAtDepthMax);
    return i >= 0
        || ((i = -i - 1) < polyVerticesHash.length && (polyVerticesHash[i] >> twiceDeltaDepth) == hash) 
        || (i > 0 && (polyVerticesHash[--i] >> twiceDeltaDepth) == hash);
  }

  private static final CooXYZ[] buildArrayOfCooXYZ(double[][] vertices) {
    final CooXYZ[] polyVertices = new CooXYZ[vertices.length];
    for (int i = 0; i < vertices.length; i++) {
      final double[] lonlat = vertices[i];
      polyVertices[i] = new CooXYZ(lonlat[LON_INDEX], lonlat[LAT_INDEX]);
    }
    return polyVertices;
  }

  private final long[] buildArrayOfPolyVerticesHash(final CooXYZ[] polyVertices) {
    final long[] polyVerticesHash = new long[polyVertices.length];
    for (int i = 0; i < polyVertices.length; i++) {
      final CooXYZ polyVertex = polyVertices[i];
      polyVerticesHash[i] = this.hcDepthMax.hash(polyVertex.lon(), polyVertex.lat());
    }
    return polyVerticesHash;
  }

}
