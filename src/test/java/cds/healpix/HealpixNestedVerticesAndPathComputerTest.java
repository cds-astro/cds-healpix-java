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

import java.io.File;
import java.io.IOException;
/*import java.nio.file.Files;
import java.nio.file.StandardOpenOption;*/
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.HealpixNested;
import cds.healpix.VerticesAndPathComputer;
import cds.healpix.CompassPoint.Cardinal;

public class HealpixNestedVerticesAndPathComputerTest {

  @Test
  public void centerAndVerticesDepth0Test() throws IOException {
    centerAndVertices(0);
  }
  
  @Test
  public void centerAndVerticesDepth1Test() throws IOException {
    centerAndVertices(1);
  }
  
  @Test
  public void centerAndVerticesDepth2Test() throws IOException {
    centerAndVertices(2);
  }
  
  @Test
  public void centerAndVerticesDepth3Test() throws IOException {
    centerAndVertices(3);
  }
  
  @Test
  public void centerAndVerticesDepth4Test() throws IOException {
    centerAndVertices(4);
  }
  
  @Test
  public void centerAndVerticesDepth5Test() throws IOException {
    centerAndVertices(5);
  }
  
  public void centerAndVertices(int depth) throws IOException {
    final HealpixNested hn = Healpix.getNested(depth);
    final VerticesAndPathComputer vpc = hn.newVerticesAndPathComputer();
    final List<String> rows = new ArrayList<String>((int) hn.nHash * (1 + 4));
    rows.add("h,lonDeg,latDeg,isCenter");
    for (int i = 0; i < hn.nHash; i++) {
      final double[] lonlat = vpc.center(i);
      rows.add(lonlat2txt(i, lonlat, true));
      final EnumMap<Cardinal, double[]> lonlats = vpc.vertices(i, EnumSet.allOf(Cardinal.class));
      for (final double[] ll : lonlats.values()) {
        rows.add(lonlat2txt(i, ll, false));
      }
    }
    final File f = new File("./target/test-results/centersAndVertices.depth" + depth + ".csv");
    /*Files.write(f.toPath(), rows, StandardOpenOption.CREATE, StandardOpenOption.WRITE,
        StandardOpenOption.TRUNCATE_EXISTING);*/
  }
  
  private String lonlat2txt(final long hash, final double[] coos, boolean isCenter) {
    return String.format(Locale.US, "%d,%.10f,%.10f,%s",
        hash, Math.toDegrees(coos[0]), Math.toDegrees(coos[1]), isCenter);
  }
  
}
