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

import java.util.EnumMap;

import cds.healpix.CompassPoint.Cardinal;
import cds.healpix.common.math.HackersDelight;

import static cds.healpix.CompassPoint.Cardinal.E;
import static cds.healpix.CompassPoint.Cardinal.N;
import static cds.healpix.CompassPoint.Cardinal.S;
import static cds.healpix.CompassPoint.Cardinal.W;
import static cds.healpix.Healpix.TRANSITION_Z;
import static cds.healpix.HealpixNestedBMOC.buildValue;
import static cds.healpix.HealpixUnprojector.ONE_OVER_SQRT6;
import static cds.healpix.Projection.LAT_INDEX;
import static cds.healpix.Projection.LON_INDEX;
import static cds.healpix.common.math.FastMath.acos;
import static cds.healpix.common.math.FastMath.asin;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.ONE_OVER_2PI;
import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.TWO_PI;
import static cds.healpix.common.math.Math.abs;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;


/**
 * Remarque:
 *  - at the poles, cone parts which are not in the cell containing the cone center cannot be part
 *    of a cell if the cell do not have a vertex in the cone!
 *  - At NPC and SPC, if baseCell != baseCellCenter, then take the symmetric of NE or NW (or SE or SW).
 *  
 *  Si coneCenterCell = pole (i.e. (d0h < 4 || 7 < d0h) && x=nide-1 && y=nside-1),
 *  regarder deltaLon max:
 *   deltaLonMax = arcsin(sinCenterLat / cosRadius)   (does not works well for small angles, near the poles)
 *   
 * 
 * @author F.-X. Pineau
 *
 */
final class NestedSmallCell implements HealpixNestedFixedRadiusConeComputer {

//  private static final Logger LOG = LoggerFactory.getLogger("fr.unistra.cds.healpix.cone"); // NestedLargeCell.class.getName()

  private final int nIterMax;
  private final double relativePrecision;

  private final int startingDepth;
  private final int deeperDepth;
  private final int deltaDepthMax;
  private final HealpixNested[] hns;
  private final VerticesAndPathComputer[] vpComputers;

  private final int startingNside;
  private final int startingNsideTime3;
  private final int startingNsideTime4;
  
  private final double rRad;
  private final double twoSineOfHalfConeRadius;
  private final double squareOfsinOfHalfR;

  private final HashComputer hComputerStartingDepth;
  private final HashComputer hComputerDepthMax;
  
  private final NeighbourSelector neigSelector;
  
  private final ConeOrdinalHashComputer cohc;
  
  private final AngularDistanceComputer angDistComputer;

  private final double[] center = new double[2];
  private final EnumMap<Cardinal, double[]> vertices = new EnumMap<Cardinal, double[]>(Cardinal.class);

  private final FlatHashList neigList;

  private final long[] specialHashs = new long[4]; 
  private int iSpecialHash;

  private int baseCellHash;             // Hash value of the base cell (i.e. the depth 0 cell)
  private int iInBaseCell, jInBaseCell; // South-East / South-West coordinates inside a base cell

  private int smallestCornerRingIndex;
  
  //Because we do not want HealpixNestedLonLatComputer to publicly implement SettableHashPartsstartingNsideTime3
  private final SettableHashParts hashPartsProxy = new SettableHashParts() {
    @Override public int baseCellHash() { return baseCellHash; }
    @Override public int iInBaseCell()  { return iInBaseCell; }
    @Override public int jInBaseCell()  { return jInBaseCell; }
    @Override public void setBaseCellHash(int baseCelHash) { baseCellHash = baseCelHash; }
    @Override public void setIInBaseCell(int iInBaseCel)   { iInBaseCell = iInBaseCel; }
    @Override public void setJInBaseCell(int jInBaseCel)   { jInBaseCell = jInBaseCel; }
  };

  private static interface AdditionalCheck {
    boolean isOk(final int deltaDepth, final long hash, final double[] deltasLon,
        final double coneCenterLonRad);
  }
  private static final AdditionalCheck ALWAYS_OK = new AdditionalCheck() {
    @Override
    public boolean isOk(final int deltaDepth, final long hash, final double[] deltasLon,
        final double coneCenterLonRad) {
      return true;
    }
  };
  private final AdditionalCheck centerInCone = new AdditionalCheck() {
    @Override
    public boolean isOk(final int deltaDepth, final long hash, final double[] deltasLon,
        final double coneCenterLonRad) {
      return NestedSmallCell.this.containsCenter(deltaDepth, hash, deltasLon, coneCenterLonRad);
    }
  };
  
  public NestedSmallCell(final int startingDepth, final int deeperDepth,
      final double radiusRad, final ConeOrdinalHashComputer cohc) {
    this(startingDepth, deeperDepth, radiusRad, cohc, radiusRad * 1e-4, 10);
  }

  protected NestedSmallCell(final int startingDepth, final int deeperDepth,
      final double radiusRad, final ConeOrdinalHashComputer cohc, final double eps, final int nIterMax) {
    assert radiusRad > 0;
    this.startingDepth = startingDepth;
    this.deeperDepth = deeperDepth;
    this.deltaDepthMax = this.deeperDepth - this.startingDepth;
    this.hns = new HealpixNested[this.deltaDepthMax + 1];
    this.vpComputers = new VerticesAndPathComputer[this.deltaDepthMax + 1];
    for (int i = 0 ; i <= this.deltaDepthMax; i++) {
      final HealpixNested hn = Healpix.getNested(this.startingDepth + i);
      this.hns[i] = hn;
      this.vpComputers[i] = hn.newVerticesAndPathComputer();
    }
    this.rRad = radiusRad;
    this.relativePrecision = this.rRad * eps;
    this.nIterMax = nIterMax;
    this.angDistComputer = AngularDistanceComputer.getComputer(this.rRad);
    final double sinHalfR = angDistComputer.sin(0.5 * this.rRad); 
    this.twoSineOfHalfConeRadius = 2 * sinHalfR;
    this.squareOfsinOfHalfR = sinHalfR * sinHalfR;

    this.hComputerDepthMax = this.hns[this.deltaDepthMax].newHashComputer();
    
    final HealpixNested hn0 = this.hns[0];
    this.hComputerStartingDepth = hn0.newHashComputer();
    this.neigSelector = hn0.newNeighbourSelector();
    this.startingNside = Healpix.nside(this.startingDepth);
    this.startingNsideTime3 = hn0.nsideTime(3);
    this.startingNsideTime4 = hn0.nsideTime(4);
    
    this.cohc = cohc;
    
    this.neigList = new FlatHashList(-1, 9); // We don't care, internal usage only
    
    initVerticesMap();
  }

  private final void initVerticesMap() {
    for (final Cardinal c : Cardinal.values()) {
      vertices.put(c, new double[2]);
    }
  }

  protected static final int nElemMax(final int nRings, final int nRingOtherSideOfThePole) {
    //   nRing * 2      (smallest cells at the extremity of each ring)
    // + nRing / 2 * 2  (cells of depth deppestDepth - 1 a the extremity of each ring)
    // + nRing / 4 * 2
    // + ...
    // + nRing / (2^DeltaDepth) * 2
    // = nRing * (2 + 1 + 1/2 + 1/4 + ...)
    // and the serie 1/2 + 1/4 + ... = 1
    // ~= nRing * (2 + 1 + 1) = 4  * nRing
    // Constant should be 4, we use 6 to be conservative
    return 6 * (nRings + nRingOtherSideOfThePole);
  }
  
  
  @Override
  public final double getRadius() {
    return this.rRad;
  }

  @Override
  public HealpixNestedFixedRadiusConeComputer newComputer() {
    return new NestedSmallCell(this.startingDepth, this.deeperDepth, this.rRad, this.cohc);
  }
  
  @Override
  public final HealpixNestedBMOC overlappingCells(double coneCenterLonRad, double coneCenterLatRad) {
    return overlapping(coneCenterLonRad, coneCenterLatRad, ALWAYS_OK);
  }

  @Override
  public HealpixNestedBMOC overlappingCenters(double coneCenterLonRad, double coneCenterLatRad) {
    return overlapping(coneCenterLonRad, coneCenterLatRad, this.centerInCone);
  }
  
  public final HealpixNestedBMOC overlapping(double coneCenterLonRad, double coneCenterLatRad, final AdditionalCheck check) {
    // Pre-compute constants
    final double cosConeCenterLat = cos(coneCenterLatRad);
    final double sinConeCenterLat = sin(coneCenterLatRad);
    coneCenterLonRad = normalizeLon(coneCenterLonRad);
    assert 0 <= coneCenterLonRad && coneCenterLonRad <= TWO_PI;

    // Compute hash of the cell containing the cone center
    final long centerHash = this.hComputerStartingDepth.hash(coneCenterLonRad, coneCenterLatRad);
    assert -HALF_PI <= coneCenterLatRad && coneCenterLatRad <= HALF_PI;

    // Count the number of rings in the other sie of a pole (if the cone contain a pole),
    // the smallest possible ring index...
    int nRingsOtherSideOfPole = 0;
    this.smallestCornerRingIndex = -1;  
    double latMax = coneCenterLatRad + this.rRad;
    if (latMax > HALF_PI) {
      nRingsOtherSideOfPole = ringIndex(this.deltaDepthMax, this.hComputerDepthMax.hash(coneCenterLonRad + PI, HALF_PI - latMax)) + 1;
      latMax = HALF_PI;
    } else {
      this.smallestCornerRingIndex // min = -1 = North polevpComputers
      = Math.max(this.smallestCornerRingIndex,
          ringIndex(this.deltaDepthMax, this.hComputerDepthMax.hash(coneCenterLonRad, latMax)) - 1);
    }
    // ... and the largest possible ring index
    final int nIsolatRings = Healpix.nIsolatitudeRings(this.deeperDepth);
    int largestCornerRingIndex = nIsolatRings;
    double latMin = coneCenterLatRad - this.rRad;
    if (latMin < -HALF_PI) {
      nRingsOtherSideOfPole = largestCornerRingIndex - ringIndex(this.deltaDepthMax, this.hComputerDepthMax.hash(coneCenterLonRad + PI, -HALF_PI - latMin));
      latMin = -HALF_PI;
    } else {
      largestCornerRingIndex = Math.min(largestCornerRingIndex, ringIndex(this.deltaDepthMax, this.hComputerDepthMax.hash(coneCenterLonRad, latMin)) + 1); 
    }
// LOG.info("Ring indexes. lo: {}; hi: {}; max: {}.", smallestCornerRingIndex, largestCornerRingIndex, Healpix.nIsolatitudeRings(this.deeperDepth));
    // Deduce the number of rings
    final int nRings = largestCornerRingIndex - this.smallestCornerRingIndex + 1;
    // Compute deltaLon for each ring
    final double[] deltasLon = new double[nRings];
    for (int i = 0; i < nRings; i++) {
      final double ringLat = latOf(this.deltaDepthMax, this.smallestCornerRingIndex + i);
// LOG.info("ringLat {}; deltasLon {}.", ringLat);
      deltasLon[i] = deltaLon(ringLat - coneCenterLatRad, cos(ringLat), cosConeCenterLat); // Double.NaN if not in the cone
    }
    if (smallestCornerRingIndex == -1) {
      deltasLon[0] = PI;
    }
    if (largestCornerRingIndex == nIsolatRings) {
      deltasLon[nRings - 1] = PI;
    }
 // LOG.info("nRings {}; deltasLon {}.", nRings, Arrays.toString(deltasLon));
    // Now compute special points hashs
    this.iSpecialHash = this.cohc.computeOrdinalHash(coneCenterLonRad, coneCenterLatRad, this.rRad,
        this.hComputerDepthMax, this.angDistComputer, this.relativePrecision, this.nIterMax,
        centerHash, this.vpComputers[0], this.vertices,
        cosConeCenterLat, sinConeCenterLat, this.twoSineOfHalfConeRadius, this.squareOfsinOfHalfR,
        this.specialHashs);
    // Now start the work on all cells recursively :)
    this.neigSelector.neighbours(centerHash, this.neigList);
    this.neigList.put(centerHash);
    this.neigList.sortByHashAsc();
    final long[] mocElems = new long[nElemMax(nRings, nRingsOtherSideOfPole)];
    int mocSize = 0;
    for (int i = 0; i < this.neigList.size(); i++) {
// LOG.info("Test: {}", this.neigList.get(i));
      mocSize = buildMocRecursively(mocElems, mocSize, 0, this.neigList.get(i), deltasLon, coneCenterLonRad, check);
    }
    return HealpixNestedBMOC.createUnsafe(this.deeperDepth, mocElems, mocSize);
  }
  
  private final int buildMocRecursively(final long[] moc, int mocLength,
      int deltaDepth, long hash, final double[] deltasLon,
      final double coneCenterLonRad, final AdditionalCheck check) {
    final int nV = nVerticesInLonRanges(deltaDepth, hash, deltasLon, coneCenterLonRad);
// System.out.println("dd: " + deltaDepth + "; hash: " + hash + "; nV: " + nV);
    if (nV == 4) {
      moc[mocLength++] = buildValue(this.startingDepth + deltaDepth, hash, true, this.deeperDepth);
    } else if (nV > 0 || isInSpecialHash(this.deltaDepthMax - deltaDepth, hash)) {
      if (deltaDepth == this.deltaDepthMax) {
        assert this.startingDepth + deltaDepth == this.deeperDepth;
// System.out.println(" - 1: dd: " + deltaDepth + "; hash: " + hash + "; nV: " + nV);
        if (check.isOk(deltaDepth, hash, deltasLon, coneCenterLonRad)) {
// System.out.println(" - 2: dd: " + deltaDepth + "; hash: " + hash + "; nV: " + nV);
          moc[mocLength++] = buildValue(this.deeperDepth, hash, false, this.deeperDepth);
        }
      } else {
        hash <<= 2;
        deltaDepth++;
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth,   hash, deltasLon, coneCenterLonRad, check);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, deltasLon, coneCenterLonRad, check);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, deltasLon, coneCenterLonRad, check);
        mocLength = buildMocRecursively(moc, mocLength, deltaDepth, ++hash, deltasLon, coneCenterLonRad, check);
      }
    }
    return mocLength;
  }
  
  private final int nVerticesInLonRanges(final int deltaDepth, final long hash, final double[] deltasLon,
      final double coneCenterLonRad) {
    final int dd = this.deltaDepthMax - deltaDepth; // assert this.startingDepth + detaDepth; 
    final int nside = 1 << dd;
    final int ringIndexMax = ringIndex(this.deltaDepthMax, hash << (dd << 1)) + 1 - this.smallestCornerRingIndex;
    final int ringIndexCenter = ringIndexMax - nside;
    final int ringIndexMin = ringIndexCenter - nside ;
    this.vpComputers[deltaDepth].vertices(hash, this.vertices);
    int n = 0;
    if (0 <= ringIndexMin && ringIndexMin < deltasLon.length) {
      n += oneIfIsInCone(coneCenterLonRad, this.vertices.get(N), deltasLon[ringIndexMin]);
    }
    if (0 <= ringIndexCenter && ringIndexCenter < deltasLon.length) {
      n += oneIfIsInCone(coneCenterLonRad, this.vertices.get(E), deltasLon[ringIndexCenter]);
      n += oneIfIsInCone(coneCenterLonRad, this.vertices.get(W), deltasLon[ringIndexCenter]);
    }
    if (0 <= ringIndexMax && ringIndexMax < deltasLon.length) {
       n += oneIfIsInCone(coneCenterLonRad, this.vertices.get(S), deltasLon[ringIndexMax]);
    }
    return n;  
  }
  
  private final boolean containsCenter(final int deltaDepth, final long hash, final double[] deltasLon,
      final double coneCenterLonRad) {
    final int dd = this.deltaDepthMax - deltaDepth; // assert this.startingDepth + detaDepth; 
    final int nside = 1 << dd;
    final int ringIndexMax = ringIndex(this.deltaDepthMax, hash << (dd << 1)) + 1 - this.smallestCornerRingIndex;
    final int ringIndexCenter = ringIndexMax - nside;
    this.vpComputers[deltaDepth].center(hash, center);
    return isInCone(coneCenterLonRad, center, deltasLon[ringIndexCenter]);
  }
  
  private static final int oneIfIsInCone(final double coneCenterLonRad, final double[] lonLat,
      final double deltaLonMax) {
    return isInCone(coneCenterLonRad, lonLat, deltaLonMax) ? 1 : 0;
  }
  
  private static final boolean isInCone(final double coneCenterLonRad, final double[] lonLat,
      final double deltaLonMax) {
// System.out.println("coneCenterLonRad: " + coneCenterLonRad + "; lonLat[LON_INDEX]: " + lonLat[LON_INDEX] + "; deltaLonMax: " + deltaLonMax);
// System.out.println("res: " + (computeDeltaLon(coneCenterLonRad, lonLat[LON_INDEX]) <= deltaLonMax));
    return computeDeltaLon(coneCenterLonRad, lonLat[LON_INDEX]) <= deltaLonMax;
  }
  
  private final boolean isInSpecialHash(final int deltaDepth, final long hash) {
    //assert deltaDepth <= this.deltaDepthMax;
    assert this.iSpecialHash == 4;
    final int twiceDeltaDepth = deltaDepth << 1;
    // Unrolled for loop
    return (this.specialHashs[0] >>> twiceDeltaDepth) == hash
        || (this.specialHashs[1] >>> twiceDeltaDepth) == hash
        || (this.specialHashs[2] >>> twiceDeltaDepth) == hash
        || (this.specialHashs[3] >>> twiceDeltaDepth) == hash;
  }

  private final double squareOfSinOfHalfOf(final double angleRad) {
    final double sin = angDistComputer.sin(0.5 * angleRad);
    return sin * sin;
  }

  /**
   * Return a value between 0 and 1, or NaN if deltaLatRad is larger than the radius of the cone.
   * We recall that any comparison (<, >, <=, >=, ==) other than (!=)  with NaN returns false 
   * @param deltaLatRad
   * @param cosLat
   * @param cosLatCenter
   * @return
   */
  private final double squareOfSinOfHalfDeltaLon(final double deltaLatRad,
      final double cosLat, final double cosLatCenter) {
    assert (cosLat >= 0 && cosLatCenter >= 0) || Double.isNaN(cosLat) : cosLat + " " + cosLatCenter;
    final double n = this.squareOfsinOfHalfR - squareOfSinOfHalfOf(deltaLatRad);
    final double d = cosLat * cosLatCenter;
    if (n < 0) { //  || !Double.isFinite(n) || !Double.isFinite(d)
      return Double.NaN;
    } else if (n >= d) {
      return 1; 
    } else {
      return n / d;  // return NaN if n or d is NaN
    }
  }

  private final double deltaLon(final double deltaLatRad,
      final double cosLat, final double cosLatCenter) {
// LOG.info("deltaLatRad {}; cosLat: {}; cosLatCenter: {}.", deltaLatRad, cosLat, cosLatCenter);
    return 2 * Math.asin(Math.sqrt(
        squareOfSinOfHalfDeltaLon(deltaLatRad, cosLat, cosLatCenter)
        ));
  }

  protected static final double normalizeLon(final double lonRad) {
    if (lonRad < 0 || TWO_PI < lonRad) {
      final int k = HackersDelight.floorInt(lonRad * ONE_OVER_2PI);
      return lonRad - k * TWO_PI;
    }
    return lonRad;
  }

  private static final double computeDeltaLon(final double lonARad, final double lonBRad) {
    final double res = abs(lonARad - lonBRad);
    return res <= PI ? res : TWO_PI - res;
  }

 private final int ringIndex(final int deltaDepth, final long hash) {
    final HealpixNested hn = this.hns[deltaDepth];
    hn.decodeRegularHash(hash, this.hashPartsProxy);
    return ringIndex(hn, this.baseCellHash, this.iInBaseCell, this.jInBaseCell);
  }

  private static final int ringIndex(final HealpixNested hn,
      long baseCellHash, int iInBasePixel, int jInBasePixel) {
    final int h = iInBasePixel + jInBasePixel;
    final long jBasePixel = hn.dividedBy4Quotient(baseCellHash);
    return (int) (hn.nsideTime(jBasePixel + 2) - (h + 2));
  }

  final double latOf(final int deltaDepth, final int iRing) {
    if (iRing < 0 || iRing > (this.startingNsideTime4 << deltaDepth)) {
      return Double.NaN;
    }
    // C.f. HealpixUnprojector and method decodeRegularRing()
    // final double y = 2 - ((iRing + 1) / (double) hn.nside);
    final int nside = this.startingNside << deltaDepth;
    double lat = ((iRing + 1) / (double) nside);
    if (isInNorthPolarCap(iRing, nside - 1)) { // North polar cap
      lat *= ONE_OVER_SQRT6;
      lat = 2 * acos(lat) - HALF_PI;
    } else if (isInSouthPolarCap(iRing, this.startingNsideTime3 << deltaDepth)) { // South polar cap
      lat = 4 - lat;
      lat *= ONE_OVER_SQRT6;
      lat = -2 * acos(lat) + HALF_PI;
    } else { // Equatorial region
      lat = asin((2 - lat) * TRANSITION_Z);
    }
    return lat;
  }

  private static final boolean isInNorthPolarCap(final int iRing, final int nsideMinus1) {
    return iRing < nsideMinus1;
  }

  private static final boolean isInSouthPolarCap(final int iRing, final int nsideTime3) {
    return iRing >= nsideTime3;
  }

}
