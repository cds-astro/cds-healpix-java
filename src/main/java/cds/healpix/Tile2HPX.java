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

import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

import cds.healpix.fillingcurve.FillingCurve2D;
import cds.healpix.fillingcurve.FillingCurve2DType;
import cds.healpix.fillingcurve.ZOrderCurve2DImpls;

import static cds.healpix.common.math.Math.PI_OVER_FOUR;

/**
 * 
 * @author F.-X. Pineau
 *
 */
public final class Tile2HPX {
  
  public enum WCSFrame {
    EQU ("RA--", "DEC-"),
    GAL ("GLON", "GLAT");

    private final String lonLabel;
    private final String latLabel;

    private WCSFrame(final String lonLabel, final String latLabel) {
      this.lonLabel = lonLabel;
      this.latLabel = latLabel;
    }
  }

  public final FillingCurve2D fc = ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_INT;
  
  // Input parameters
  private final int order;
  private final int inNside;
  private final WCSFrame frame;
  // Precomputed quantities
  private final int nsideTile;
  private final int nsidePix;
  private final int twiceDepth;
  private final long xyMask;
  

  public Tile2HPX(final int tileOrder, final int inTileNside, final WCSFrame wcsFrame) {
    FillingCurve2DType.Z_ORDER_LUPT.get(tileOrder);
    this.order = tileOrder;
    this.inNside = inTileNside;
    this.frame = wcsFrame;
    // Derivated quantities
    this.nsideTile = 1 << this.order;
    this.nsidePix = this.nsideTile * this.inNside;
    this.twiceDepth = tileOrder << 1;
    this.xyMask = (1L << this.twiceDepth) - 1;
  }

  /**
   * Returns the coordinates of the center of the cell defined by the given hash in the HEALPix
   * projection plane.
   * @param hash value defining the cell we want the coordinate of the center
   * @param xyResult object used to store the result
   */
  public void center(long hash, final double[] xyResult) {
    // Pull apart the hash elements
    final int d0h = (int) (hash >> this.twiceDepth);
    hash = this.fc.hash2ij(hash & this.xyMask);
    final int iInD0h = this.fc.ij2i(hash);
    final int jInD0h = this.fc.ij2j(hash);
    // Compute coordinates from the center of the base pixel with  x-axis = W-->E, y-axis = S-->N
    final int lInD0h = iInD0h - jInD0h;
    final int hInD0h = iInD0h + jInD0h - (this.nsideTile - 1);
    // Compute coordinates of the base pixel in the projection plane
    final int d0hBy4Quotient = d0h >> 2;
    final int d0hMod4 = d0h - (d0hBy4Quotient << 2);
    final int hD0h = 1 - d0hBy4Quotient;
    int lD0h = d0hMod4 << 1;
    if ((hD0h == 0 && (lD0h == 6 || (lD0h == 4 && lInD0h > 0))) // case equatorial region
        || (hD0h != 0 && ++lD0h > 3)) { // case polar caps regions
      lD0h -= 8;
    }
    // Finalize computation
    xyResult[0] = PI_OVER_FOUR * (lD0h + lInD0h / (double) nsideTile);
    xyResult[1] = PI_OVER_FOUR * (hD0h + hInD0h / (double) nsideTile);
  }

 
  public final Map<String, String> toFitsHeader(final long tileIpix) throws Exception {
    final double[] xy = new double[2];
    this.center(tileIpix, xy);

    final double centreX = Math.toDegrees(xy[0]);
    final double centreY = Math.toDegrees(xy[1]);

    final double scale = 45d / nsidePix;

    // cooAvantDeproj = [CD][cooPixel - CRPIX]
    // cooPixel - (inNside + 1) / 2.0 => put origin at the center of the tile
    // [CD][T] = [DeltaX, DeltaY] => T = [CD]^-1  [DeltaX, DeltaY]

    final double crPix1 = +((inNside + 1) / 2.0) - 0.5 * (-centreX / scale + centreY / scale);
    final double crPix2 = +((inNside + 1) / 2.0) - 0.5 * (-centreX / scale - centreY / scale);

    final Map<String, String> output = new LinkedHashMap<String, String>();
    output.put("NAXIS  ", "                    2 / number of data axes");
    output.put("NAXIS1 ", String.format(Locale.US, "%21d / length of data axis 1", inNside));
    output.put("NAXIS2 ", String.format(Locale.US, "%21d / length of data axis 1", inNside));

    output.put("CRPIX1 ", String.format(Locale.US, "%21.1f / Coordinate reference pixel", crPix1));
    output.put("CRPIX2 ", String.format(Locale.US, "%21.1f / Coordinate reference pixel", crPix2));

    // Solution using CDs
    output.put("CD1_1  ", String.format(Locale.US, "%21.13E / Transformation matrix (rot + scale)", -scale));
    output.put("CD1_2  ", String.format(Locale.US, "%21.13E / Transformation matrix (rot + scale)", -scale));
    output.put("CD2_1  ", String.format(Locale.US, "%21.13E / Transformation matrix (rot + scale)", +scale));
    output.put("CD2_2  ", String.format(Locale.US, "%21.13E / Transformation matrix (rot + scale)", -scale));

    /*// Solution using PC and CDELT
    output.put("PC1_1  ", String.format(Locale.US, "%21.13f / Transformation matrix (rot + scale)", -0.5 * Math.sqrt(2)));
    output.put("PC1_2  ", String.format(Locale.US, "%21.13f / Transformation matrix (rot + scale)", -0.5 * Math.sqrt(2)));
    output.put("PC2_1  ", String.format(Locale.US, "%21.13f / Transformation matrix (rot + scale)", +0.5 * Math.sqrt(2)));
    output.put("PC2_2  ", String.format(Locale.US, "%21.13f / Transformation matrix (rot + scale)", -0.5 * Math.sqrt(2)));
    output.put("CDELT1 ", String.format(Locale.US, "%21.13f / [deg] Coordinate increment", scale / (0.5 * Math.sqrt(2))));
    output.put("CDELT2 ", String.format(Locale.US, "%21.13f / [deg] Coordinate increment", scale / (0.5 * Math.sqrt(2))));
     */

    output.put("CTYPE1 ", " '" + frame.lonLabel + "-HPX'           / Longitude in an HPX projection");
    output.put("CTYPE2 ", " '" + frame.latLabel + "-HPX'           /  Latitude in an HPX projection");

    output.put("CRVAL1 ", "                   0. / [deg] Longitude at the reference point");
    output.put("CRVAL2 ", "                   0. / [deg]  Latitude at the reference point");
    output.put("PV2_1  ", "                   4  / HPX H parameter (longitude)");
    output.put("PV2_2  ", "                   3  / HPX K parameter  (latitude)");

    return output;
  }


  public static final Map<String, String> toFitsHeader(final int tileOrder,
      final long tileIpix, final int nPixelX, final WCSFrame frame) throws Exception {
    return new Tile2HPX(tileOrder, nPixelX, frame).toFitsHeader(tileIpix);
  }

  public static void main(final String[] args) throws Exception {
    int order = Integer.parseInt(args[0]);
    long ipix = Long.parseLong(args[1]);
    int nPixelsX = Integer.parseInt(args[2]);

    // final Map<String, String> map = toFitsHeader(8, 300373L/*289450*/, 512, WCSFrame.EQU);
    final Map<String, String> map =  toFitsHeader(order, ipix, nPixelsX,
        args.length > 3 ? WCSFrame.valueOf(args[3].toUpperCase()): WCSFrame.EQU);
    for (final Map.Entry<String, String> e : map.entrySet()) {
      System.out.println(e.getKey() + " =" + e.getValue());
    }
  }

}
