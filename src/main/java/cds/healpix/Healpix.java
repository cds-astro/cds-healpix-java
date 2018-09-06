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

import cds.healpix.common.math.FastMath;
import cds.healpix.fillingcurve.FillingCurve2DType;

import static cds.healpix.common.math.FastMath.acos;
import static cds.healpix.common.math.FastMath.sinQ;
import static cds.healpix.common.math.HackersDelight.BUT_SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.isPowOf2Fast;
import static cds.healpix.common.math.HackersDelight.toBits;
import static cds.healpix.common.math.Math.FOUR_OVER_PI;
import static cds.healpix.common.math.Math.HALF_PI;
import static cds.healpix.common.math.Math.PI;
import static cds.healpix.common.math.Math.PI_OVER_FOUR;
import static cds.healpix.common.math.Math.SQRT6;
import static cds.healpix.common.math.Math.asin;
import static cds.healpix.HealpixUnprojector.EPS_POLE;
import static cds.healpix.HealpixUnprojector.ONE_OVER_SQRT6;;


/**
 * Utility class containing HEALPix constants and operations which are independent of the depth.
 *
 * We replace the term "norder" used in HEALPix documents by the more generic "depth" usually used
 * in tree-datastructures. We consider each base-resolution cell as been part of the root so that
 * their depth equals 0.
 *
 * @author F.-X. Pineau
 *
 */
public final class Healpix implements Projection {

  /** 29, largest possible depth we can store on a signed positive long
   * (4 bits for base cells + 2 bits per depth + 2 remaining bits (1 use in the unique notation). */
  public static final int DEPTH_MAX = 29;
  private static final int NSIDE_MAX = nside(DEPTH_MAX);

  /** Gorski2005, Section 4: "Number of base-resolution pixel layers between the N and S poles".
   * Written N_\Theta in Gorsky2005 and K in Calabretta2007. */
  private static final int N_LAT = 3;

  /** Gorski2005, Section 4: "Number of equatorial [...] base-resolution pixels".
   * Written N_\Phi in Gorsky2005 and H in Calabretta2007. */
  private static final int N_LON = 4;

  /** Gorski2005, Section 4: "Number of base-resolution pixels", i.e. of depth 0 hash values.
   * Written N_{PIX} in Gorsky2007. */
  private static final int N_DEPTH0_HASH = N_LAT * N_LON;

  /** Limit on |z|=|sin(lat)| between the equatorial region and the polar caps.
   * Equals 2/3, see Eq. (1) in Gorsky2005. */
  public static final double TRANSITION_Z = (N_LAT - 1) / (double) N_LAT; // = 2/3 Noted Z_THR ?
  static final double ONE_OVER_TRANSITION_Z = 1.5d; // = 1 / TRANSITION_Z;
  
  /** Limit on the |latitude (in radians)| between the equatorial region and the polar caps.
   * Equals asin(2/3) = 0.7297276562269663 radians ~= 41,81 degrees.
   * Written \theta_X in Calabretta2007. */
  public static final double TRANSITION_LATITUDE = asin(TRANSITION_Z);

  /**
   * For each HEALPix depth, stores the smallest distance from an edge of a cell to the opposite
   * edge of the same cell. If the radius of a cone is smaller than this distance, we know that
   * it will overlap maximum 9 pixels (the pixel containing the center of the cone plus
   * the 8 neighbours).
   * In practice, this distance if the distance between the point of coordinate
   * (0, {@link Healpix#TRANSITION_LATITUDE}) and it nearest point on the Northeast edge of the
   * cell of base hash 0 and coordinates in the base hash (x=0, y=nside-1).
   * IMPORTANT REMARK:
   * - this value is larger than the smallest center to vertex distance
   * - this value x2 is larger than the smallest diagonal (NS or EW)
   * - this value x2 is larger than the smallest edge
   * - BUT there is no case in which the value is larger than the four center-to-vertex distance
   * - BUT there is no case in which the value x2 is larger than both diagonals
   * =&gt; this radius is smaller than the smaller circumcircle radius (=&gt; no cone having the smaller
   * -edge-to-opposite-edge-radius radius can contains the 4 vertices of a cell (but 3 is ok)
   * vertices 
   */
  private static double[] SMALLER_EDGE2OPEDGE_DIST = new double[] {
      0.8410686705685088,    // depth = 0
      0.37723631722170053,   // depth = 1
      0.18256386461918295,   // depth = 2
      0.09000432499034523,   // depth = 3
      0.04470553761855741,   // depth = 4
      0.02228115704023076,   // depth = 5
      0.011122977211214961,  // depth = 6
      0.005557125022105058,  // depth = 7
      0.0027774761500209185, // depth = 8
      0.0013884670480328143, // depth = 9
      6.941658374603201E-4,  // depth = 10
      3.4706600585087755E-4, // depth = 11
      1.7352877579970442E-4, // depth = 12
      8.676333125510362E-5,  // depth = 13
      4.338140148342286E-5,  // depth = 14
      2.1690634707822447E-5, // depth = 15
      1.084530084565172E-5,  // depth = 16
      5.422646295795749E-6,  // depth = 17
      2.711322116099695E-6,  // depth = 18
      1.3556608000873442E-6, // depth = 19
      6.778303355805395E-7,  // depth = 20
      3.389151516386149E-7,  // depth = 21
      1.69457571754776E-7,   // depth = 22
      8.472878485272006E-8,  // depth = 23
      4.236439215502565E-8,  // depth = 24
      2.1182195982014308E-8, // depth = 25
      1.0591097960375205E-8, // depth = 26
      5.295548939447981E-9,  // depth = 27
      2.647774429917369E-9,  // depth = 28
      1.3238871881399636E-9};// depth = 29

  
  /**
   * Latitude, in the equatorial region, for which the distance from the cell center to its four
   * vertices is equal on the sky (i.e. the shape of the cell on the sky is similar to a square):
   * dX = dY = 1 / nside (center to vertex distance)
   * X = 4/pi * lon     =&gt; dX = 4/pi dlon
   * Y = 3/2 * sin(lat) =&gt; dY = 3/2 * cos(lat) dlat
   * dlon * cos(lat) = dlat (same distance on the sky)
   * =&gt; cos^2(lat) = 2/3 * 4/pi
   * =&gt; lat = arccos(sqrt(2/3 * 4/pi)) = 22.88050802005781976715 deg = 0.39934019947897773410 rad
   */
  public static final double LAT_OF_SQUARE_CELL = 0.39934019947897773410;
  static final double COS_LAT_OF_SQUARE_CELL = 0.92131773192356127804;
  //public static final double TWICE_COS_LAT_OF_SQUARE_CELL = 2 * COS_LAT_OF_SQUARE_CELL;
  
  /** Constant = (1 - COS_LAT_OF_SQUARE_CELL) /  LAT_OF_SQUARE_CELL^2 . */
  private static final double EQUAT_PARABOLA_COEFF = 0.4933905296766163;
  /** Constant = pi / 4 * 1 / nside.
   * It is the distance between the Center of a cell and its East (or West) vertex on the equator. */
  private static final double[] DIST_CW_EQUATOR = new double[]{
      0.7853981633974483,    // depth = 0
      0.39269908169872414,   // depth = 1
      0.19634954084936207,   // depth = 2
      0.09817477042468103,   // depth = 3
      0.04908738521234052,   // depth = 4
      0.02454369260617026,   // depth = 5
      0.01227184630308513,   // depth = 6
      0.006135923151542565,  // depth = 7
      0.0030679615757712823, // depth = 8
      0.0015339807878856412, // depth = 9
      7.669903939428206E-4,  // depth = 10
      3.834951969714103E-4,  // depth = 11
      1.9174759848570515E-4, // depth = 12
      9.587379924285257E-5,  // depth = 13
      4.7936899621426287E-5, // depth = 14
      2.3968449810713143E-5, // depth = 15
      1.1984224905356572E-5, // depth = 16
      5.992112452678286E-6,  // depth = 17
      2.996056226339143E-6,  // depth = 18
      1.4980281131695715E-6, // depth = 19
      7.490140565847857E-7,  // depth = 20
      3.7450702829239286E-7, // depth = 21
      1.8725351414619643E-7, // depth = 22
      9.362675707309822E-8,  // depth = 23
      4.681337853654911E-8,  // depth = 24
      2.3406689268274554E-8, // depth = 25
      1.1703344634137277E-8, // depth = 26
      5.8516723170686385E-9, // depth = 27
      2.9258361585343192E-9, // depth = 28
      1.4629180792671596E-9  // depth = 29
  };
  /** Constant = asin(1 - (1 - 1 / nside)^2 / 3) - TRANSITION_LATITUDE.
   * It is the distance between the transition latitude ring (Y = 1) and the previous 
   * ring (Y = 1 + 1 / nside).
   * So, it is the distance between the Center and the North vertex of a cell having its center
   * on the transition latitude, and the its north pole on the vertical. */
  private static final double[] DIST_TRANSIT_LAT_TO_PREV_RING = new double[]{
      0.8410686705679302,    // depth = 0
      0.4299308082455824,    // depth = 1
      0.21870018201290964,   // depth = 2
      0.11049492567137786,   // depth = 3
      0.05556369564623209,   // depth = 4
      0.027864942299874795,  // depth = 5
      0.013953769172063923,  // depth = 6
      0.006982275959440676,  // depth = 7
      0.0034924942664104064, // depth = 8
      0.0017465872656792225, // depth = 9
      8.733787988972619E-4,  // depth = 10
      4.367107076018728E-4,  // depth = 11
      2.1836068292035993E-4, // depth = 12
      1.0918167400042478E-4, // depth = 13
      5.459117016770598E-5,  // depth = 14
      2.7295668379889726E-5, // depth = 15
      1.3647855014453647E-5, // depth = 16
      6.823932713562186E-6,  // depth = 17
      3.411967658295545E-6,  // depth = 18
      1.7059841546096521E-6, // depth = 19
      8.529921586841738E-7,  // depth = 20
      4.2649609977019054E-7, // depth = 21
      2.1324805488109888E-7, // depth = 22
      1.0662402882832822E-7, // depth = 23
      5.331201469171987E-8,  // depth = 24
      2.6656007512393387E-8, // depth = 25
      1.3328003811707845E-8, // depth = 26
      6.664001905853922E-9,  // depth = 27
      3.33200089741581E-9,   // depth = 28
      1.666000448707905E-9   // depth = 29
  };
  
  /**
   * 
   */
  private static final double[] POLAR_LINEAR_COEFF = new double[]{
      Double.NaN,            // depth = 0
      0.1312192099038766,    // depth = 1
      0.060492837392452155,  // depth = 2
      0.029001500919391094,  // depth = 3
      0.01419382425997025,   // depth = 4
      0.007020691494838553,  // depth = 5
      0.0034913519041639258, // depth = 6
      0.0017409349692292947, // depth = 7
      8.692831611333698E-4,  // depth = 8
      4.343456141662374E-4,  // depth = 9
      2.1709882974034925E-4, // depth = 10
      1.0853092231337608E-4, // depth = 11
      5.426083823966203E-5,  // depth = 12
      2.7129263418298822E-5, // depth = 13
      1.356434278725635E-5,  // depth = 14
      6.782099163568464E-6,  // depth = 15
      3.391031524325326E-6,  // depth = 16
      1.6955112477746438E-6, // depth = 17
      8.477544953072421E-7,  // depth = 18
      4.238769654916027E-7,  // depth = 19
      2.1193841220245034E-7, // depth = 20
      1.0596918844808161E-7, // depth = 21
      5.298458983383655E-8,  // depth = 22
      2.649229379051603E-8,  // depth = 23
      1.324614662808308E-8,  // depth = 24
      6.62307322705194E-9,   // depth = 25
      3.311536588893444E-9,  // depth = 26
      1.6557682911737157E-9, // depth = 27
      8.278841563091061E-10, // depth = 28
      4.1394207794998953E-10 // depth = 29
  };
  /**
   * xa = LAT_OF_SQUARE_CELL;
   * xb = TRANSITION_LATITUDE;
   * ya = CST_B[d] * COS_LAT_OF_SQUARE_CELL;
   * yb = CST_C[d];
   * EQUAT_LINEAR_COEFF[d] = (yb - ya) / (xb - xa)
   */
  private static double[] EQUAT_LINEAR_COEFF = new double[]{
      0.35554441795671393,   // depth = 0
      0.20621297683954326,   // depth = 1
      0.11441072474544801,   // depth = 2
      0.06067048988147291,   // depth = 3
      0.03129240237539165,   // depth = 4
      0.01589770733619387,   // depth = 5
      0.008013317441437975,  // depth = 6
      0.004022977056353704,  // depth = 7
      0.0020155936679524915, // depth = 8
      0.0010088263297003488, // depth = 9
      5.046709411546653E-4,  // depth = 10
      2.523999650150198E-4,  // depth = 11
      1.2621611241593917E-4, // depth = 12
      6.311208947320893E-5,  // depth = 13
      3.155705315107901E-5,  // depth = 14
      1.5778778691759518E-5, // depth = 15
      7.889452376446927E-6,  // depth = 16
      3.944741946495324E-6,  // depth = 17
      1.9723749126056044E-6, // depth = 18
      9.861884413943153E-7,  // depth = 19
      4.930944670120405E-7,  // depth = 20
      2.4654729533676847E-7, // depth = 21
      1.232736627900346E-7,  // depth = 22
      6.163683559547574E-8,  // depth = 23
      3.081841863782956E-8,  // depth = 24
      1.5409209822969793E-8, // depth = 25
      7.704605079503234E-9,  // depth = 26
      3.852302539751617E-9,  // depth = 27
      1.926151101857471E-9,  // depth = 28
      9.630755509287354E-10  // depth = 29
  };
  
  /**
   * For a given position on the unit sphere and a given depth, returns an upper limit and the
   * distance between the center of a cell and its farthest vertex.
   * @param lonRad longitude in radians.
   * @param latRad latitude in radians.
   * @param depth HEALPix depth. 
   * @return 2/3 * 1/nside * cos(LAT_OF_SQUARE_CELL)
   */
  public static double getLargestCenterToCellVertexDistance(
      final double lonRad, double latRad, final int depth) {
    latRad = Math.abs(latRad);
    if (depth == 0) {
      return latRad < LAT_OF_SQUARE_CELL ? PI_OVER_FOUR : HALF_PI - TRANSITION_LATITUDE;
    }
    // dX = dY = 1 / nside (center to vertex distance)
    // Cylindrical projection:
    //   X = 4/pi * lon     => dX = 4/pi dlon
    //   Y = 3/2 * sin(lat) => dY = 3/2 * cos(lat) dlat
    // Collignon projection:
    //   Y = sqrt( 3 (1-sin(lat)) ) => lat = asin(1 - Y^2 / 3)
    // At LAT_OF_SQUARE_CELL (same distance on the sky):
    //   dlon * cos(lat) = dlat, i.e. d(CE) = d(CW) = 4/pi * 1/nside * cos(LAT_OF_SQUARE_CELL) ~= d(CN) ~= d(CS)
    // At lat = 0:
    //   d(CE) = d(CW) = dlon = pi/4 * 1/nside
    if (latRad < LAT_OF_SQUARE_CELL) {
      // For lat < LAT_OF_SQUARE_CELL, dlon > dlat
      // Parabola approximation from lat = 0 to lat = LAT_OF_SQUARE_CELL, such that
      //        a * LAT_OF_SQUARE_CELL^2 + b = 4/pi * 1/nside * cos(LAT_OF_SQUARE_CELL)
      //   and  a * 0^2 + b = pi/4 * 1/nside
      //   => b = pi/4 * 1/nside
      //   => a = 4/pi * 1/nside (1 - cos(LAT_OF_SQUARE_CELL)) / LAT_OF_SQUARE_CELL^2
      //        = b * [ (1 - cos(LAT_OF_SQUARE_CELL)) / LAT_OF_SQUARE_CELL^2 ]
      return DIST_CW_EQUATOR[depth] * (1 - EQUAT_PARABOLA_COEFF * (latRad * latRad)); // Parabola approximation
    } else if (latRad < TRANSITION_LATITUDE) {
      // Linear aprox between A and B with
      //   A: lat = LAT_OF_SQUARE_CELL, d(CE) = d(CW) =  4/pi * 1/nside * cos(LAT_OF_SQUARE_CELL)
      //   B: lat = TRANSITION_LATITUDE, d(CN) = asin(1 - ((1 - 1/nside)^2 / 3)) - TRANSITION_LATITUDE 
      //      i.e. North vertex at the vertical with Y = 2 - (1 + 1/nside) 
      // Then simple linear approx
      return DIST_CW_EQUATOR[depth] * COS_LAT_OF_SQUARE_CELL + (latRad - LAT_OF_SQUARE_CELL) * EQUAT_LINEAR_COEFF[depth];
    } else {
      // For lat >= TRANSITION_LATITUDE
      // Linear approx from d(CN) at TRANSITION_LATITUDE, border of a base cell, and d(CN) at
      // TRANSITION_LATITUDE and center of a base cell.
      //   At transition latitude, border
      //     Cell Center,  C: X = 1/nside, Y = 2 - 1             => (lon = pi/4 * 1 / nside, lat = TRANSITION_LATITUDE)
      //     North vertex, N: X = 0,       Y = 2 - (1 + 1/nside) => (lon = 0               , lat = asin(1 - ((1 - 1/nside)^2 / 3))
      //     We note the distance CN
      //   At transition latitude, center
      //     d(CN) = asin(1 - (1 - 1 / nside)^2 / 3) - TRANSITION_LATITUDE
      // Then make a simple linear approximation:
      // CST_D = (CN - d(CN)) / ( PI/4  * (1 - 1d / nside) )
      final double lonMod = Math.abs(PI_OVER_FOUR - lonRad % HALF_PI);
      return POLAR_LINEAR_COEFF[depth] * lonMod + DIST_TRANSIT_LAT_TO_PREV_RING[depth];
    }
  }
  
  /**
   * Default implementation of the z-order curve. Set to the look-up table implementation because it
   * is the one showing the best performances in simple tests.<br>
   * For an application requiring a lot of CPU cache, you may try the Z_ORDER_XOR implementation.
   */
  public static final FillingCurve2DType DEFAULT_FCTYPE = FillingCurve2DType.Z_ORDER_LUPT;
  
  private static final HealpixNested[] HEALPIX_NESTED = new HealpixNested[30];

  private static final HealpixNestedFast[] HEALPIX_NESTED_FAST = new HealpixNestedFast[30];

  private static final EnumMap<FillingCurve2DType, HealpixNested[]> FCTYPE_HEALPIX_NESTED = 
      new EnumMap<FillingCurve2DType, HealpixNested[]>(FillingCurve2DType.class);
      
  private static final EnumMap<FillingCurve2DType, HealpixNestedFast[]> FCTYPE_HEALPIX_NESTED_FAST = 
          new EnumMap<FillingCurve2DType, HealpixNestedFast[]>(FillingCurve2DType.class);
  
  static {
    for (final FillingCurve2DType fct : FillingCurve2DType.values()) {
      if (fct == DEFAULT_FCTYPE) {
        FCTYPE_HEALPIX_NESTED.put(fct, HEALPIX_NESTED);
        FCTYPE_HEALPIX_NESTED_FAST.put(fct, HEALPIX_NESTED_FAST);
      } else {
        FCTYPE_HEALPIX_NESTED.put(fct, new HealpixNested[30]);
        FCTYPE_HEALPIX_NESTED_FAST.put(fct, new HealpixNestedFast[30]);
      }
    }
  }
  
  /** Unique instance of this class (needed to implement the {@link Projection} interface). */
  public static final Healpix UI = new Healpix();
  
  /** Prevents from instantiation out of the internal unique instance. */
  private Healpix() { }

  /**
   * Lazy instantiation of unique instances of {@link HealpixNested} for each required depth.
   * Uses internally a double-checked lock.
   * @param depth the depth of the wanted {@link HealpixNested} instance
   * @return a {@link HealpixNested} instance at the given depth
   */
  public static HealpixNested getNested(final int depth) {
    HealpixNested instance = HEALPIX_NESTED[depth];
    if (instance == null) {
      synchronized(HEALPIX_NESTED) {
        instance = HEALPIX_NESTED[depth];
        if (instance == null) {
          instance = new HealpixNested(depth, DEFAULT_FCTYPE);
          HEALPIX_NESTED[depth] = instance;
        }
      }
    }
    return instance;
  }
  
  // getRing ??
  // getNestedWest ??
  // getPeanoHilbert ??
  
  /**
   * Same as {@link #getNested(int)} except that we have here a given instance per depth and 
   * z-order curve implementation {@link #DEFAULT_FCTYPE}.
   * @param depth the depth of the wanted {@link HealpixNested} instance
   * @param fillingCurveType the wanted z-order curve implementation
   * @return a {@link HealpixNested} instance at the given depth, with the given z-order curve
   *         implementation
   */
  public static HealpixNested getNested(final int depth, FillingCurve2DType fillingCurveType) {
    final HealpixNested[] hn = FCTYPE_HEALPIX_NESTED.get(fillingCurveType);
    HealpixNested instance = hn[depth];
    if (instance == null) {
      synchronized(HEALPIX_NESTED) {
        instance = hn[depth];
        if (instance == null) {
          instance = new HealpixNested(depth, fillingCurveType);
          hn[depth] = instance;
        }
      }
    }
    return instance;
  }

  /**
   * Same as {@link #getNested(int)} for {@link HealpixNestedFast}, i.e. faster but less readable
   * version of {@link HealpixNested}.
   * @param depth the depth of the wanted {@link HealpixNestedFast} instance
   * @return a {@link HealpixNestedFast} instance at the given depth
   */
  public static HealpixNestedFast getNestedFast(final int depth) {
    HealpixNestedFast instance = HEALPIX_NESTED_FAST[depth];
    if (instance == null) {
      synchronized(HEALPIX_NESTED_FAST) {
        instance = HEALPIX_NESTED_FAST[depth];
        if (instance == null) {
          instance = new HealpixNestedFast(depth, DEFAULT_FCTYPE);
          HEALPIX_NESTED_FAST[depth] = instance;
        }
      }
    }
    return instance;
  }
  
  /**
   * Same as {@link #getNested(int, FillingCurve2DType)} for {@link HealpixNestedFast},
   * i.e. faster but less readable version of {@link HealpixNested}.
   * @param depth he depth of the wanted {@link HealpixNestedFast} instance
   * @param fillingCurveType the wanted z-order curve implementation
   * @return a {@link HealpixNestedFast} instance at the given depth, with the given z-order curve
   *         implementation
   */
  public static HealpixNestedFast getNestedFast(final int depth, FillingCurve2DType fillingCurveType) {
    final HealpixNestedFast[] hn = FCTYPE_HEALPIX_NESTED_FAST.get(fillingCurveType);
    HealpixNestedFast instance = hn[depth];
    if (instance == null) {
      synchronized(HEALPIX_NESTED_FAST) {
        instance = hn[depth];
        if (instance == null) {
          instance = new HealpixNestedFast(depth, DEFAULT_FCTYPE);
          hn[depth] = instance;
        }
      }
    }
    return instance;
  }
  
  /**
   * Returns the number of cells along both axis of a base-resolution cell at the given depth.
   * @param depth (or order) number of subdivision of a base-resolution cell, from 0 to {@link #DEPTH_MAX}. 
   * @return 2^depth.
   */
  public static int nside(final int depth) {
    checkDepth(depth);
    return 1 << depth; // i.e. 2^depth
  }

  // Used to mutliply a double by nside/2
  static long halfNside4IEEEdouble(int depth) { 
    return depth == 0 ? -(1L << 52) : (depth - 1L) << 52;
  }
  // Used to mutliply a double by nside
  static long nside4IEEEdouble(int depth) { 
    return depth << 52;
  }
  
  static void checkDepth(final int depth) {
    if (depth < 0 || depth > DEPTH_MAX) {
      throw new IllegalArgumentException("Expected depth in [0, " + DEPTH_MAX + "]. Actual: " + depth);
    }
  }

  static void checkLatitude(final double latRad) {
    if (latRad < -HALF_PI || HALF_PI < latRad) {
      throw new IllegalArgumentException("Wrong latitude. Expected: in [-pi/2, pi/2]. Actual: " + latRad);
    }
  }

  /**
   * Returns the number of subdivision of a base-resolution pixel at the given nside.
   * @param nside number of pixels along both axis of a base-resolution pixel, from 1 to 2 power 
   *              {@link #DEPTH_MAX}.
   * @return log2(nside)
   */
  public static int depth(final int nside) {
    if (!isPowOf2Fast(nside) || nside < 1 || nside > NSIDE_MAX) {
      throw new IllegalArgumentException("Nside must be a power of 2 in [1-2^" + DEPTH_MAX + "]");
    }
    return Integer.numberOfTrailingZeros(nside);
  }

  static long nsideSquare(final int depth) {
    checkDepth(depth);
    return 1L << (depth << 1);
  }
  
  /**
   * Returns the number of cells (so the number of distinct hash value) the unit sphere is divided
   * in at the given depth.
   * @param depth (or order) number of subdivision of a base-resolution cell, from 0 to {@link #DEPTH_MAX}.
   * @return 12 * nside^2
   */
  public static long nHash(final int depth) {
    checkDepth(depth);
    return ((long) N_DEPTH0_HASH) << (depth << 1); // i.e. 12 * nside^2
  }  
  /**
   * Returns the number of isolatitude rings at the given depth, i.e. the number of small circles
   * parallel to the equator containing HEALPix cell centers
   * @param depth (or order) number of subdivision of a base-resolution cell, from 0 to {@link #DEPTH_MAX}.
   * @return 4 * nside - 1
   */
  public static int nIsolatitudeRings(final int depth) {
    // 2 * (2 * nside - 1) + 1 , i.e. 2 * number of isolatitude rings per polar cap cell
    // + 1 in the equator
    checkDepth(depth);
    assert (4 << depth) == 4 * nside(depth);
    return (4 << depth) - 1;    
  }

  /**
   * Returns {@code true} if the given latitude is in the North polar cap
   * @param latRad latitude, in radians
   * @return {@code true} if {@code latRad >} {@link #TRANSITION_LATITUDE}
   */
  public static boolean isLatInNorthPolarCap(final double latRad) {
    return latRad > TRANSITION_LATITUDE;
  }
  
  /**
   * Returns {@code true} if the given latitude is in the South polar cap
   * @param latRad latitude, in radians
   * @return {@code true} if {@code latRad < -}{@link #TRANSITION_LATITUDE}
   */
  public static boolean isLatInSouthPolarCap(final double latRad) {
    return latRad < -TRANSITION_LATITUDE;
  }

  static boolean isInEquatorialRegion(final double absLatRad) {
    return absLatRad <= TRANSITION_LATITUDE;
  }

  /**
   * Returns the the smallest depth (in [0, 29]) at which a shape having the given largest distance
   * from its center to a border overlaps a maximum of 9 cells (the cell containing the center of
   * the shape plus the 8 neighbouring cells).<br>
   * Info: internally, unrolled binary search loop on 30 pre-computed values (one by depth).
   * @param distMaxInRad largest possible distance, in radians, between the center and the border
   *        of a shape.
   * @return -1 if the given distance is very large (&gt; ~48deg), else returns the smallest depth 
   * (in [0, 29]) at which a shape having the given largest distance from its center to a border
   * overlaps a maximum of 9 cells (the cell containing the center of the shape plus the 8 
   * neighbouring cells).
   */
  public static int getBestStartingDepth(final double distMaxInRad) {
    // Unrolled binary search loop
    if (distMaxInRad > SMALLER_EDGE2OPEDGE_DIST[0]) {
      return -1;
    } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[29]) {
      return 29;
    } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[15]) {
      // in [15, 28]
      if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[22]) {
        // in [22, 28]
        if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[25]) {
          // in [25, 28]
          if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[27]) {
            // in [27, 28]
            if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[28]) {
              return 28;
            } else {
              return 27;
            } // in [25, 26]
          } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[26]) {
            return 26;
          } else {
            return 25;
          } // in [22, 24]
        } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[24]) {
          return 24;
        } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[23]) {
          return 23;
        } else {
          return 22;
        } // in [15, 21]
      } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[18]) {
        // in [18, 21]
        if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[20]) {
          // in [20, 21]
          if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[21]) {
            return 21;
          } else {
            return 20;
          } // in [18, 19]
        } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[19]) {
          return 19;
        } else {
          return 18;
        } // in [15, 17]
      } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[17]) {
        return 17;
      } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[16]) {
        return 16;
      } else {
        return 15;
      } // in [0, 14]
    } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[7]) {
      // in [7, 14]
      if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[11]) {
        // in [11, 14]
        if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[13]) {
          // in [13, 14]
          if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[14]) {
            return 14;
          } else {
            return 13;
          } // in [11, 12]
        } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[12]) {
          return 12;
        } else {
          return 11;
        } // in [7, 10]
      } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[9]) {
        // in [9, 10]
        if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[10]) {
          return 10;
        } else {
          return 9;
        } // in [7, 8]
      } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[8]) {
        return 8;
      } else {
        return 7;
      } // in [0, 6]
    } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[3]) {
      // in [3, 6]
      if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[5]) {
        // in [5, 6]
        if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[6]) {
          return 6;
        } else {
          return 5;
        } // in [3, 4]
      } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[4]) {
        return 4;
      } else {
        return 3;
      }  // in [0, 2]
    } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[2]) {
      return 2;
    } else if (distMaxInRad < SMALLER_EDGE2OPEDGE_DIST[1]) {
      return 1;
    } else {
      return 0;
    }
  }

  
  static double sqrtOfThreeTimeOneMinusSinOf(double a) {
    // Eq. 19 in Calabretta (2007)
    // sqrt(3 * (1 - sin(x)))
    // = sqrt(3 * (1 + cos(x + pi/2)))
    // = sqrt(3 * 2 * (1 + cos(2 * (x/2 + pi/4))) / 2 )
    // = sqrt(6) * cos(x/2 + pi/4)
    assert 0 <= a && a <= HALF_PI;
    return SQRT6 * FastMath.cosQ(0.5 * a + PI_OVER_FOUR);
  }

  /**
   * Project the given spherical coordinates into the Euclidean plane, following the HPX WCS
   * projection (see Calabretta2007).
   * Remark: this method is thread safe
   * @param lonRad longitude in radians, support positive and negative reasonably large values with a
   *        naive approach (no Cody-Waite nor Payne Hanek range reduction).
   * @param latRad latitude, must be in [-pi/2, pi/2] radians
   * @return x in [-PI, PI] radians, y in [-PI/2, PI/2] radians.
   */
  @Override
  public double[] project(double lonRad, double latRad) {
    final double[] res = new double[2];
    project(lonRad, latRad, res);
    return res;
  }

  @Override
  public void project(double lonRad, double latRad, double[] resultXY) {
    // See class HealpixProjector for a cleaner code!!
    checkLatitude(latRad);
    // compute absolute value of lon and keep its sign
    long signLon = toBits(lonRad);
    final double absLon = fromBits(signLon & BUT_SIGN_BIT_MASK_L);
    signLon &= SIGN_BIT_MASK_L;
    // compute absolute value of lat and keep its sign
    long signLat = toBits(latRad);
    final double absLat = fromBits(signLat & BUT_SIGN_BIT_MASK_L);
    signLat &= SIGN_BIT_MASK_L;
    // scale absLon in [-1, 1] and keep range offset (in n*PI/4)
    double x = FOUR_OVER_PI * absLon, y;
    int xOffset = ((int) x) | 1;
    x -= xOffset;
    xOffset &= 7; // reduce range offset to max 7;
    if (isInEquatorialRegion(absLat)) { // Cylindrical Equal Area projection
      y = sinQ(absLat, 0);
      y *= ONE_OVER_TRANSITION_Z;
    } else { // Collignon projection
      y = sqrtOfThreeTimeOneMinusSinOf(absLat);
      x *= y; 
      y = (2 - y);
    }
    // We do this to obtain values in [-pi, pi] instead of possibly [0,2pi]
    xOffset -= (xOffset & 4) << 1; // <=> xOffset = xOffset < 4 ? xOffset : xOffset - 8;
    x += xOffset;
    x = fromBits(signLon | toBits(x)); // apply lon sign to x
    y = fromBits(signLat | toBits(y)); // apply lat sign to y
    resultXY[X_INDEX] = x * PI_OVER_FOUR;
    resultXY[Y_INDEX] = y * PI_OVER_FOUR;
  }

  @Override
  public double[] unproject(double x, double y) {
    final double[] res = new double[2];
    unproject(x, y, res);
    return res;
  }

  @Override
  public void unproject(double x, double y, double[] resultLonLat) {
    // See class HealpixUnprojector for a cleaner code!!
    checkProjectionBounds(x, y); // return NaN, NaN instead of throwing an error?
    // compute |x| and keep the x sign bit
    long signX = toBits(x);
    final double absX = fromBits(signX & BUT_SIGN_BIT_MASK_L) * FOUR_OVER_PI;
    signX &= SIGN_BIT_MASK_L;
    // compute |y| and keep the y sign bit
    long signY = toBits(y);
    final double absY = fromBits(signY & BUT_SIGN_BIT_MASK_L) * FOUR_OVER_PI;
    signY &= SIGN_BIT_MASK_L;
    // shift |x| in [-1, +1] and keep the offset
    // int xFloor = (int) absX;
    int xOffset = ((int) absX) | 1;
    double lon = absX - xOffset,  lat;
    if (isInPlaneEquatorialRegion(absY)) {
      lat = asin(absY * TRANSITION_Z);
    } else {
      lat = 2 - absY;
      if (isNotNearFromPole(lat)) { // Rare, so few risks of branch miss-prediction
        lon /= lat;
        lon = dealWithNumericalApproxInEdge(lon);
      } // in case of pole, lon = lat = 0 (we avoid NaN due to division by lat=0)
      lat *= ONE_OVER_SQRT6;
      lat = 2 * acos(lat) - HALF_PI;
    }
    // shift lon by offset and scale to radians
    lon += xOffset;
    // apply x sign to lon;
    lon = fromBits(signX | toBits(lon));
    lon += (signX >>> 60); // <=> +8 to have only positive values
    lon *= PI_OVER_FOUR;
    // apply y sign to lat;
    lat = fromBits(signY | toBits(lat));
    // store results in the given array
    resultLonLat[LON_INDEX] = lon;
    resultLonLat[LAT_INDEX] = lat;
  } 

  private static void checkProjectionBounds(final double x, final double y) {
    if (x < -PI || x > PI) {
      throw new IllegalArgumentException("x value \"" + x + "\" must be in [-pi, pi]");
    }
    if (y < -HALF_PI || y > HALF_PI) {
      throw new IllegalArgumentException("y value \"" + x + "\" must be in [-pi / 2, pi / 2]");
    }
  }
  private static boolean isInPlaneEquatorialRegion(double absY) {
    return absY <= 1;
  }
  private static boolean isNotNearFromPole(final double sqrtOfThreeTimeOneMinusSinOf) {
    // In case of pole: x = y = 0
    return sqrtOfThreeTimeOneMinusSinOf > EPS_POLE;
  }
  private double dealWithNumericalApproxInEdge(double lon) { // Rare, so few risks of branch miss-prediction
    return lon > 1 ? 1 : lon < -1 ? -1 : lon;
  }
  
  /**
   * Create the unique representation of the given hash at the given depth.
   * The unique representation encode the depth together with the hash value such that each possible
   * (deph, hash) pair is unique.<br>
   * To do so, the unique representation uses a sentinel bit to code the depth.
   * The sentinel bit is the (1 +4 + 2*depth)^th most significant bit.<br>
   * The encoding in the case of the nested scheme is thus<br>
   * <b> 0...0sbbbb112233...</b><br>
   * With<br>
   * <ul>
   *   <li>0...0: unused bits</li>
   *   <li>s: sentinel bit</li>
   *   <li>bbbb: the 4 bits coding the base cell</li>
   *   <li>11: the 2 bits coding depth 1</li>
   *   <li>22: the 2 bits coding depth 2</li> 
   *   <li>33: the 2 bits coding depth 3</li> 
   *   <li>...</li> 
   * </ul>
   * @param depth the depth of the wanted unique hash
   * @param hash the hash we want the unique representation
   * @return the unique representation of the given hash at the given depth.
   */
  public static long uniq(final int depth, final long hash) {
    return (16L << (depth << 1)) | hash;
  }
  
  /**
   * Extract the depth from the unique representation of the hash.
   * See {@link #uniq(int, long)} to know more about the uniq encoding.
   * @param uniqedHash the uniq nested hash value we want to extract the depth
   * @return the depth from the unique representation of the hash.
   */
  public static int uniq2depth(final long uniqedHash) {
    return (60 - Long.numberOfLeadingZeros(uniqedHash)) >>> 1;
  }

  /**
   * Extract the hash from the the unique representation of the hash.
   * See {@link #uniq(int, long)} to know more about the uniq encoding.
   * @param uniqedHash the uniq nested hash value we want to extract the hash value
   * @return the hash from the the unique representation of the hash.
   */
  public static long uniq2hash(final long uniqedHash) {
    return  uniqedHash & ~Long.highestOneBit(uniqedHash);
  }
  
  /**
   * Faster version of {@link #uniq2hash(long)} in case we know (or have already extracted) the depth.
   * @param uniqdHash the uniq hash value we want to extract the hash value
   * @param depth the known depth encoded in the given uniq hash
   * @return the hash from the the unique representation of the hash, knowing the depth.
   */
  public static long uniq2hash(final long uniqdHash, final int depth) {
    return uniqdHash & ~(16L << (depth << 1));
  }

  /*private static final void genCST_A() {
    double v = (1 - COS_LAT_OF_SQUARE_CELL) /  (LAT_OF_SQUARE_CELL * LAT_OF_SQUARE_CELL);
    System.out.println("double EQUAT_PARABOLA_COEFF = " + v + ";");
  }
  
  private static final void genCST_B() {
    System.out.println("double[] DIST_CW_EQUATOR = new double[]{");
    for (int d = 0; d <= Healpix.DEPTH_MAX; d++) {
      double val = PI_OVER_FOUR / Healpix.nside(d);
      System.out.print(val);
      if (d != Healpix.DEPTH_MAX) {
        System.out.println(", // depth = " + d);
      }
    }
    System.out.println("};");
  }
  
  private static final void genCST_C() {
    AngularDistanceComputer adc = AngularDistanceComputer.getComputer(1.0);
    System.out.println("double[] DIST_TRANSIT_LAT_TO_PREV_RING = new double[]{");
    for (int d = 0; d <= Healpix.DEPTH_MAX; d++) {
      final int nside = Healpix.nside(d);
      double val = 1 - 1d / nside;
      double lat = asin(1 - (val * val) / 3d);
      val = lat - TRANSITION_LATITUDE;
      System.out.print(val);
      if (d != Healpix.DEPTH_MAX) {
        System.out.println(", // depth = " + d);
      }
    }
    System.out.println("};");
  }
  
  private static final void genCST_D() {
    System.out.println("double[] POLAR_LINEAR_COEFF = new double[]{");
    AngularDistanceComputer adc = AngularDistanceComputer.getComputer(1.0);
    for (int d = 0; d <= Healpix.DEPTH_MAX; d++) {
      // Largest distance (lonMod = +-pi/4)
      // => Center: X= 1/nside, Y = 1           => (lon = pi/4 * 1 / nside, lat = TRANSITION_LATITUDE)
      // => North vertex: X = 0, Y =1 + 1/nside => (lon = 0, lat = asin[1 - (1 - 1/nside)^2/3]
      // => Y = 1 - sqrt(3*(1-sin(lat))) => lat = asin[1 - (1 - Y)^2/3]
      // Collignon:
      // Y = 1 - Z
      // Z = sqrt(3 * (1 - z))  in [  1, 0]
      // z = sin(lat)           in [2/3, 1]
      // => lat = arsin(z) = asin( 1 - (1 - Y)^2 / 3 )
      // => z = 1 - Z^2 / 3 = 1 - (1 - Y)^2 / 3
      // => z = 1 - (1 - 1/nside)^2 / 3
      // Y = (1 / nside)
      final int nside = Healpix.nside(d);
      final double lonA = 0.25 * PI / nside;
      final double latA = TRANSITION_LATITUDE;
      final double lonB = 0;
      double latB = 1 - 1d / nside;
      latB = asin (1  - (latB * latB) / 3d);
      // System.out.println("latB: " + latB + "; (lonB - lonA): " + (lonA) + "; (latB - latA): " + (latB - latA));
      double val = adc.haversineDistInRad(lonA, latB - latA, Math.cos(latA), Math.cos(latB));
      // val = largest angular distance
      // Compute the strait line coeff
      val = (val - DIST_TRANSIT_LAT_TO_PREV_RING[d]) / (0.25 * PI  * (1 - 1d / nside));
      System.out.print(val);
      if (d != Healpix.DEPTH_MAX) {
        System.out.println(", // depth = " + d);
      }
    }
    System.out.println("};");
  }
  
  private static final void genCST_E() {
    System.out.println("double[] EQUAT_LINEAR_COEFF = new double[]{");
    for (int d = 0; d <= Healpix.DEPTH_MAX; d++) {
      final double xa = LAT_OF_SQUARE_CELL;
      final double xb = TRANSITION_LATITUDE;
      final double ya = DIST_CW_EQUATOR[d] * COS_LAT_OF_SQUARE_CELL;
      final double yb = DIST_TRANSIT_LAT_TO_PREV_RING[d];
      double val = ((yb - ya) / (xb - xa)) ;
      System.out.print(val);
      if (d != Healpix.DEPTH_MAX) {
        System.out.println(", // depth = " + d);
      }
    }
    System.out.println("};");
  }
  
  public static void main(final String[] args) {
    // genSMALLER_EDGE2OPEDGE_DIST();
    genCST_A();
    genCST_B();
    genCST_C();
    genCST_D();
    genCST_E();
  }*/
  
}
