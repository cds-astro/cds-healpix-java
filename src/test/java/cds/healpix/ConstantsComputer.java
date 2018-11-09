// Copyright 2017-2018 - Universite de Strasbourg/CNRS
// The CDS HEALPix library is developped by the Centre de Donnees
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

import java.io.File;
import java.io.IOException;
// import java.nio.file.Files;
// import java.nio.file.StandardOpenOption;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.LinkedHashSet;
import java.util.Locale;
import java.util.Set;

import cds.healpix.CompassPoint.Cardinal;

import static cds.healpix.CompassPoint.Cardinal.E;
import static cds.healpix.CompassPoint.Cardinal.N;
import static cds.healpix.CompassPoint.Cardinal.S;
import static cds.healpix.CompassPoint.Cardinal.W;
import static java.lang.Math.toDegrees;


/**
 * Package protected class used only to compute constants to be used inside the HEALPix code.
 * 
 * @author F.-X. Pineau
 *
 */
final class ConstantsComputer {


  /**
   * R.W. Sinnott, Virtues of the Haversine, Sky and Telescope, V.68:2, P.158, 1984
   */
  private static final double haversineRad(double ra1Rad, double dec1Rad,
      double ra2Rad, double dec2Rad){
    final double d = Math.sin(0.5 * (dec2Rad - dec1Rad));
    final double a = Math.sin(0.5 * (ra2Rad - ra1Rad));
    final double h3 = d * d + a * a * Math.cos(dec1Rad) * Math.cos(dec2Rad);
    return 2 * Math.asin(Math.sqrt(h3));
  }
  
  private static final double haversineRad(final double[] radec1Rad, final double[] radec2Rad) {
    final double d = Math.sin(0.5 * (radec2Rad[1] - radec1Rad[1]));
    final double a = Math.sin(0.5 * (radec2Rad[0] - radec1Rad[0]));
    final double h3 = d * d + a * a * Math.cos(radec1Rad[1]) * Math.cos(radec2Rad[1]);
    return 2 * Math.asin(Math.sqrt(h3));
  }
  
  
  private static void verticesDistance(final int depth) throws IOException{
    final HealpixNested hn = Healpix.getNested(depth);
    final VerticesAndPathComputer vpc = hn.newVerticesAndPathComputer();
    final EnumSet<Cardinal> directions = EnumSet.allOf(Cardinal.class);
    
    final Set<String> rows = new LinkedHashSet<String>();
    rows.add("icell,lon,lat,dCN,dCS,dCE,dCW,dMaxCV,dMaCVdir,dNS,dEW,dDiagMax,dDiagMaxDir,dNE,dSE,dSW,dNW,dEdgeMin,dEdgeMinDir,computedMaxDCV");

    for (int i = 0; i < hn.nHash; i++) {
      final double[] center = vpc.center(i);
      final EnumMap<Cardinal, double[]> vertices = vpc.vertices(i, directions);
      // Center to vertices distance
      final double distCN = haversineRad(center, vertices.get(N));
      final double distCS = haversineRad(center, vertices.get(S));
      final double distCE = haversineRad(center, vertices.get(E));
      final double distCW = haversineRad(center, vertices.get(W));
      double maxDistCV = distCN;
      int maxDistCVdir = 0;
      if (distCS > maxDistCV) {
        maxDistCV = distCS;
        maxDistCVdir = 1;
      }
      if (distCE > maxDistCV) {
        maxDistCV = distCE;
        maxDistCVdir = 2;
      }
      if (distCW > maxDistCV) {
        maxDistCV = distCW;
        maxDistCVdir = 3;
      }
      // Diagonals
      final double distNS = haversineRad(vertices.get(N), vertices.get(S));
      final double distEW = haversineRad(vertices.get(E), vertices.get(W));
      double distDiagMax = distNS;
      int maxDistDiagdir = 0;
      if (distEW > distDiagMax) {
        distDiagMax = distEW;
        maxDistDiagdir = 1;
      }
      // Edge distance
      final double distNE = haversineRad(vertices.get(N), vertices.get(E));
      final double distSE = haversineRad(vertices.get(S), vertices.get(E));
      final double distSW = haversineRad(vertices.get(S), vertices.get(W));
      final double distNW = haversineRad(vertices.get(N), vertices.get(W));
      double distEdgeMin = distNE;
      int distEdgeMinDir = 0;
      if (distSE < distEdgeMin) {
        distEdgeMin = distSE;
        distEdgeMinDir = 1;
      }
      if (distSW < distEdgeMin) {
        distEdgeMin = distSW;
        distEdgeMinDir = 2;
      }
      if (distNW < distEdgeMin) {
        distEdgeMin = distNW;
        distEdgeMinDir = 3;
      }
      final String line = toCSV(depth, i, center, 
          distCN, distCS, distCW, distCE, maxDistCV, maxDistCVdir,
          distNS, distEW, distDiagMax, maxDistDiagdir,
          distNE, distSE, distSW, distNW, distEdgeMin, distEdgeMinDir);
      rows.add(line);
    }
    
    final File f = new File("./target/test-results/distances.depth" + depth + ".csv");
    // Files.write(f.toPath(), rows, StandardOpenOption.CREATE, StandardOpenOption.WRITE,
    //    StandardOpenOption.TRUNCATE_EXISTING);
  }
  
  private static String toCSV(int depth, int icell, double[] cellCenter, 
      double dCN, double dCS, double dCE, double dCW, double dMaxCV, int dMaxCSVdir,
      double dNS, double dEW, double dDiagMax, int dDiagMaxDir,
      double dNE, double dNW, double dSE, double dSW, double dEdgeMin, int dEdgeMinDir) {
    return String.format(Locale.US, "%d,%.10f,%.10f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f",
        icell, toDegrees(cellCenter[0]), toDegrees(cellCenter[1]),
        toArcsec(dCN), toArcsec(dCS), toArcsec(dCE), toArcsec(dCW), toArcsec(dMaxCV), dMaxCSVdir,
        toArcsec(dNS), toArcsec(dEW), toArcsec(dDiagMax), dDiagMaxDir,
        toArcsec(dNE), toArcsec(dNW), toArcsec(dSE), toArcsec(dSW), toArcsec(dEdgeMin), dEdgeMinDir,
        toArcsec(Healpix.getLargestCenterToCellVertexDistance(cellCenter[0], cellCenter[1], depth)));
  }
  private static final double toArcsec(final double angle) {
    return toDegrees(angle) * 3600;
  }
  
  public static void main(final String[] args) throws IOException {
    final int depth = 8;
    verticesDistance(depth);
    
    
  }

}
