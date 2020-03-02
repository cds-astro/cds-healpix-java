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

import java.util.Arrays;

import org.junit.Test;

import cds.healpix.Healpix;
import cds.healpix.HealpixNested;
import cds.healpix.NeighbourList;
import cds.healpix.NeighbourSelector;
import cds.healpix.CompassPoint.MainWind;

public class HealpixNestedNeighbourSelectorTest {

  @Test
  public void neighboursDepth0() {
    final HealpixNested hn = Healpix.getNested(0);
    final NeighbourSelector neigSelector = hn.newNeighbourSelector();
    //  North polar cap
    // - 0
    assert neigSelector.neighbour(0L, MainWind.N)  ==  2L;
    assert neigSelector.neighbour(0L, MainWind.NE) ==  1L;
    assert neigSelector.neighbour(0L, MainWind.NW) ==  3L;
    assert neigSelector.neighbour(0L, MainWind.S)  ==  8L;
    assert neigSelector.neighbour(0L, MainWind.SE) ==  5L;
    assert neigSelector.neighbour(0L, MainWind.SW) ==  4L;
    assert neigSelector.neighbour(0L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(0L, MainWind.W)  == -1L;
    checkEquals(neigSelector.neighbours(0L), new long[]{1, 2, 3, 4, 5, 8});
    // - 1
    assert neigSelector.neighbour(1L, MainWind.N)  ==  3L;
    assert neigSelector.neighbour(1L, MainWind.NE) ==  2L;
    assert neigSelector.neighbour(1L, MainWind.NW) ==  0L;
    assert neigSelector.neighbour(1L, MainWind.S)  ==  9L;
    assert neigSelector.neighbour(1L, MainWind.SE) ==  6L;
    assert neigSelector.neighbour(1L, MainWind.SW) ==  5L;
    assert neigSelector.neighbour(1L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(1L, MainWind.W)  == -1L;
    checkEquals(neigSelector.neighbours(1L), new long[]{0, 2, 3, 5, 6, 9});
    // - 2
    assert neigSelector.neighbour(2L, MainWind.N)  ==  0L;
    assert neigSelector.neighbour(2L, MainWind.NE) ==  3L;
    assert neigSelector.neighbour(2L, MainWind.NW) ==  1L;
    assert neigSelector.neighbour(2L, MainWind.S)  == 10L;
    assert neigSelector.neighbour(2L, MainWind.SE) ==  7L;
    assert neigSelector.neighbour(2L, MainWind.SW) ==  6L;
    assert neigSelector.neighbour(2L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(2L, MainWind.W)  == -1L;
    checkEquals(neigSelector.neighbours(2L), new long[]{0, 1, 3, 6, 7, 10});
    // - 3
    assert neigSelector.neighbour(3L, MainWind.N)  ==  1L;
    assert neigSelector.neighbour(3L, MainWind.NE) ==  0L;
    assert neigSelector.neighbour(3L, MainWind.NW) ==  2L;
    assert neigSelector.neighbour(3L, MainWind.S)  == 11L;
    assert neigSelector.neighbour(3L, MainWind.SE) ==  4L;
    assert neigSelector.neighbour(3L, MainWind.SW) ==  7L;
    assert neigSelector.neighbour(3L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(3L, MainWind.W)  == -1L;
    checkEquals(neigSelector.neighbours(3L), new long[]{0, 1, 2, 4, 7, 11});
    // Equatorial region
    // - 4
    assert neigSelector.neighbour(4L, MainWind.N)  == -1L;
    assert neigSelector.neighbour(4L, MainWind.NE) ==  0L;
    assert neigSelector.neighbour(4L, MainWind.NW) ==  3L;
    assert neigSelector.neighbour(4L, MainWind.S)  == -1L;
    assert neigSelector.neighbour(4L, MainWind.SE) ==  8L;
    assert neigSelector.neighbour(4L, MainWind.SW) == 11L;
    assert neigSelector.neighbour(4L, MainWind.E)  ==  5L;
    assert neigSelector.neighbour(4L, MainWind.W)  ==  7L;
    checkEquals(neigSelector.neighbours(4L), new long[]{0, 3, 5, 7, 8, 11});
    // - 5
    assert neigSelector.neighbour(5L, MainWind.N)  == -1L;
    assert neigSelector.neighbour(5L, MainWind.NE) ==  1L;
    assert neigSelector.neighbour(5L, MainWind.NW) ==  0L;
    assert neigSelector.neighbour(5L, MainWind.S)  == -1L;
    assert neigSelector.neighbour(5L, MainWind.SE) ==  9L;
    assert neigSelector.neighbour(5L, MainWind.SW) ==  8L;
    assert neigSelector.neighbour(5L, MainWind.E)  ==  6L;
    assert neigSelector.neighbour(5L, MainWind.W)  ==  4L;
    checkEquals(neigSelector.neighbours(5L), new long[]{0, 1, 4, 6, 8, 9});
    // - 6
    assert neigSelector.neighbour(6L, MainWind.N)  == -1L;
    assert neigSelector.neighbour(6L, MainWind.NE) ==  2L;
    assert neigSelector.neighbour(6L, MainWind.NW) ==  1L;
    assert neigSelector.neighbour(6L, MainWind.S)  == -1L;
    assert neigSelector.neighbour(6L, MainWind.SE) == 10L;
    assert neigSelector.neighbour(6L, MainWind.SW) ==  9L;
    assert neigSelector.neighbour(6L, MainWind.E)  ==  7L;
    assert neigSelector.neighbour(6L, MainWind.W)  ==  5L;
    checkEquals(neigSelector.neighbours(6L), new long[]{1, 2, 5, 7, 9, 10});
    // - 7
    assert neigSelector.neighbour(7L, MainWind.N)  == -1L;
    assert neigSelector.neighbour(7L, MainWind.NE) ==  3L;
    assert neigSelector.neighbour(7L, MainWind.NW) ==  2L;
    assert neigSelector.neighbour(7L, MainWind.S)  == -1L;
    assert neigSelector.neighbour(7L, MainWind.SE) == 11L;
    assert neigSelector.neighbour(7L, MainWind.SW) == 10L;
    assert neigSelector.neighbour(7L, MainWind.E)  ==  4L;
    assert neigSelector.neighbour(7L, MainWind.W)  ==  6L;
    checkEquals(neigSelector.neighbours(7L), new long[]{2, 3, 4, 6, 10, 11});
    //  South polar cap
    // - 8
    assert neigSelector.neighbour(8L, MainWind.N)  ==  0L;
    assert neigSelector.neighbour(8L, MainWind.NE) ==  5L;
    assert neigSelector.neighbour(8L, MainWind.NW) ==  4L;
    assert neigSelector.neighbour(8L, MainWind.S)  == 10L;
    assert neigSelector.neighbour(8L, MainWind.SE) ==  9L;
    assert neigSelector.neighbour(8L, MainWind.SW) == 11L;
    assert neigSelector.neighbour(8L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(8L, MainWind.W)  == -1L;
    checkEquals(neigSelector.neighbours(8L), new long[]{0, 4, 5, 9, 10, 11});
    // - 9
    assert neigSelector.neighbour(9L, MainWind.N)  ==  1L;
    assert neigSelector.neighbour(9L, MainWind.NE) ==  6L;
    assert neigSelector.neighbour(9L, MainWind.NW) ==  5L;
    assert neigSelector.neighbour(9L, MainWind.S)  == 11L;
    assert neigSelector.neighbour(9L, MainWind.SE) == 10L;
    assert neigSelector.neighbour(9L, MainWind.SW) ==  8L;
    assert neigSelector.neighbour(9L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(9L, MainWind.W)  == -1L; 
    checkEquals(neigSelector.neighbours(9L), new long[]{1, 5, 6, 8, 10, 11});
    // - 10
    assert neigSelector.neighbour(10L, MainWind.N)  ==  2L;
    assert neigSelector.neighbour(10L, MainWind.NE) ==  7L;
    assert neigSelector.neighbour(10L, MainWind.NW) ==  6L;
    assert neigSelector.neighbour(10L, MainWind.S)  ==  8L;
    assert neigSelector.neighbour(10L, MainWind.SE) == 11L;
    assert neigSelector.neighbour(10L, MainWind.SW) ==  9L;
    assert neigSelector.neighbour(10L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(10L, MainWind.W)  == -1L;
    checkEquals(neigSelector.neighbours(10L), new long[]{2, 6, 7, 8, 9, 11});
    // - 11
    assert neigSelector.neighbour(11L, MainWind.N)  ==  3L;
    assert neigSelector.neighbour(11L, MainWind.NE) ==  4L;
    assert neigSelector.neighbour(11L, MainWind.NW) ==  7L;
    assert neigSelector.neighbour(11L, MainWind.S)  ==  9L;
    assert neigSelector.neighbour(11L, MainWind.SE) ==  8L;
    assert neigSelector.neighbour(11L, MainWind.SW) == 10L;
    assert neigSelector.neighbour(11L, MainWind.E)  == -1L;
    assert neigSelector.neighbour(11L, MainWind.W)  == -1L;
    checkEquals(neigSelector.neighbours(11L), new long[]{3, 4, 7, 8, 9, 10});
  }

  private void checkEquals(NeighbourList neig, long[] expected) {
    neig.sortByHashAsc();
    final long[] a = new long[neig.size()];
    neig.arraycopy(0, a, 0, neig.size());
    assert Arrays.equals(a, expected);
  }
  
  @Test
  public void externalEdgeTest1() {
    byte depth = 1;
    long hash = 10;
    int delta_depth = 2;
    final FlatHashList actualRes = Healpix.getNested(depth).sortedExternalEdges(hash, delta_depth);
    final long[] expectedRes = new long[]{85, 87, 93, 95, 117, 138, 139, 142, 143, 154, 176, 178, 
        184, 186, 415, 437, 439, 445, 447};
    checkEquals(actualRes, expectedRes);
  }
  
  @Test
  public void externalEdgeTest2() {
    byte depth = 1;
    long hash = 11;
    int delta_depth = 2;
    final FlatHashList actualRes = Healpix.getNested(depth).sortedExternalEdges(hash, delta_depth);
    final long[] expectedRes = new long[]{63, 95, 117, 119, 125, 127, 143, 154, 155, 158, 159, 165, 
        167, 173, 175, 239, 250, 251, 254, 255};
    checkEquals(actualRes, expectedRes);
  }

  @Test
  public void externalEdgeTest3() {
    byte depth = 0;
    long hash = 0;
    int delta_depth = 2;
    final FlatHashList actualRes = Healpix.getNested(depth).sortedExternalEdges(hash, delta_depth);
    final long[] expectedRes = new long[]{26, 27, 30, 31, 47, 53, 55, 61, 63, 69, 71, 77, 79, 90, 
        91, 94, 95, 143};
    checkEquals(actualRes, expectedRes);
  }

  private void checkEquals(FlatHashList neig, long[] expected) {
    neig.sortByHashAsc();
    final long[] a = new long[neig.size()];
    neig.arraycopy(0, a, 0, neig.size());
    assert Arrays.equals(a, expected) : "Array: " + Arrays.toString(a);
  }
  
}
