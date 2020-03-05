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

  @Test
  public void externalEdgeTest4() {
    byte depth = 3;
    long hash = 0;
    int delta_depth = 5;
    final FlatHashList actualRes = Healpix.getNested(depth).sortedExternalEdges(hash, delta_depth);
    final long[] expectedRes = new long[]{1024, 1026, 1032, 1034, 1056, 1058, 1064, 1066, 1152, 
        1154, 1160, 1162, 1184, 1186, 1192, 1194, 1536, 1538, 1544, 1546, 1568, 1570, 1576, 1578,
        1664, 1666, 1672, 1674, 1696, 1698, 1704, 1706, 2048, 2049, 2052, 2053, 2064, 2065, 2068,
        2069, 2112, 2113, 2116, 2117, 2128, 2129, 2132, 2133, 2304, 2305, 2308, 2309, 2320, 2321,
        2324, 2325, 2368, 2369, 2372, 2373, 2384, 2385, 2388, 2389, 3072, 283989, 283991, 283997,
        283999, 284021, 284023, 284029, 284031, 284117, 284119, 284125, 284127, 284149, 284151, 
        284157, 284159, 284501, 284503, 284509, 284511, 284533, 284535, 284541, 284543, 284629, 
        284631, 284637, 284639, 284661, 284663, 284669, 284671, 286037, 371370, 371371, 371374, 
        371375, 371386, 371387, 371390, 371391, 371434, 371435, 371438, 371439, 371450, 371451, 
        371454, 371455, 371626, 371627, 371630, 371631, 371642, 371643, 371646, 371647, 371690, 
        371691, 371694, 371695, 371706, 371707, 371710, 371711, 372394, 589823};
    checkEquals(actualRes, expectedRes);
  }

  
  @Test
  public void unsortedExternalEdgeTest4() {
    byte depth = 3;
    long hash = 1;
    int delta_depth = 5;
    final FlatHashList actualRes = Healpix.getNested(depth).externalEdges(hash, delta_depth);
    final long[] expectedRes = new long[]{341, 343, 349, 351, 373, 375, 381, 383, 469, 471, 477, 
        479, 501, 503, 509, 511, 853, 855, 861, 863, 885, 887, 893, 895, 981, 983, 989, 991, 1013, 
        1015, 1021, 1023, 2389, 3072, 3073, 3076, 3077, 3088, 3089, 3092, 3093, 3136, 3137, 3140, 
        3141, 3152, 3153, 3156, 3157, 3328, 3329, 3332, 3333, 3344, 3345, 3348, 3349, 3392, 3393, 
        3396, 3397, 3408, 3409, 3412, 3413, 4096, 4098, 4104, 4106, 4128, 4130, 4136, 4138, 4224, 
        4226, 4232, 4234, 4256, 4258, 4264, 4266, 4608, 4610, 4616, 4618, 4640, 4642, 4648, 4650, 
        4736, 4738, 4744, 4746, 4768, 4770, 4776, 4778, 6144, 371711, 372394, 372395, 372398, 
        372399, 372410, 372411, 372414, 372415, 372458, 372459, 372462, 372463, 372474, 372475, 
        372478, 372479, 372650, 372651, 372654, 372655, 372666, 372667, 372670, 372671, 372714, 
        372715, 372718, 372719, 372730, 372731, 372734, 372735, 375466};
    checkEquals(actualRes, expectedRes);
  }
  
  
  @Test
  public void externalEdgeBaseCellsTest() {
    byte depth = 0;
    int delta_depth = 2;
    // base cell 0
    FlatHashList actualRes = Healpix.getNested(0).sortedExternalEdges(0, delta_depth);
    long[] expectedRes = new long[]{26, 27, 30, 31, 47, 53, 55, 61, 63, 69, 71, 77, 79, 90, 91, 94, 95, 143};
    checkEquals(actualRes, expectedRes);
    // base cell 1
    actualRes = Healpix.getNested(0).sortedExternalEdges(1, delta_depth);
    expectedRes = new long[]{5, 7, 13, 15, 42, 43, 46, 47, 63, 85, 87, 93, 95, 106, 107, 110, 111, 159};
    checkEquals(actualRes, expectedRes);
    // base cell 2
    actualRes = Healpix.getNested(0).sortedExternalEdges(2, delta_depth);
    expectedRes = new long[]{15, 21, 23, 29, 31, 58, 59, 62, 63, 101, 103, 109, 111, 122, 123, 126, 127, 175};
    checkEquals(actualRes, expectedRes);
    // base cell 3
    actualRes = Healpix.getNested(0).sortedExternalEdges(3, delta_depth);
    expectedRes = new long[]{10, 11, 14, 15, 31, 37, 39, 45, 47, 74, 75, 78, 79, 117, 119, 125, 127, 191};
    checkEquals(actualRes, expectedRes);
    // base cell 4
    actualRes = Healpix.getNested(0).sortedExternalEdges(4, delta_depth);
    expectedRes = new long[]{0, 2, 8, 10, 48, 49, 52, 53, 90, 117, 138, 139, 142, 143, 181, 183, 189, 191};
    checkEquals(actualRes, expectedRes);
    // base cell 5
    actualRes = Healpix.getNested(0).sortedExternalEdges(5, delta_depth);
    expectedRes = new long[]{0, 1, 4, 5, 16, 18, 24, 26, 69, 106, 133, 135, 141, 143, 154, 155, 158, 159};
    checkEquals(actualRes, expectedRes);
    // base cell 6
    actualRes = Healpix.getNested(0).sortedExternalEdges(6, delta_depth);
    expectedRes = new long[]{16, 17, 20, 21, 32, 34, 40, 42, 85, 122, 149, 151, 157, 159, 170, 171, 174, 175};
    checkEquals(actualRes, expectedRes);
    // base cell 7
    actualRes = Healpix.getNested(0).sortedExternalEdges(7, delta_depth);
    expectedRes = new long[]{32, 33, 36, 37, 48, 50, 56, 58, 74, 101, 165, 167, 173, 175, 186, 187, 190, 191};
    checkEquals(actualRes, expectedRes);
    // base cell 8
    actualRes = Healpix.getNested(0).sortedExternalEdges(8, delta_depth);
    expectedRes = new long[]{0, 64, 65, 68, 69, 80, 82, 88, 90, 144, 146, 152, 154, 160, 176, 177, 180, 181};
    checkEquals(actualRes, expectedRes);
    // base cell 9
    actualRes = Healpix.getNested(0).sortedExternalEdges(9, delta_depth);
    expectedRes = new long[]{16, 80, 81, 84, 85, 96, 98, 104, 106, 128, 129, 132, 133, 160, 162, 168, 170, 176};
    checkEquals(actualRes, expectedRes);
    // base cell 10
    actualRes = Healpix.getNested(0).sortedExternalEdges(10, delta_depth);
    expectedRes = new long[]{32, 96, 97, 100, 101, 112, 114, 120, 122, 128, 144, 145, 148, 149, 176, 178, 184, 186};
    checkEquals(actualRes, expectedRes);
    // base cell 11
    actualRes = Healpix.getNested(0).sortedExternalEdges(11, delta_depth);
    expectedRes = new long[]{48, 64, 66, 72, 74, 112, 113, 116, 117, 128, 130, 136, 138, 144, 160, 161, 164, 165};
    checkEquals(actualRes, expectedRes);
    
    
    /*for (int h = 0; h < 12; h++) {
      actualRes = Healpix.getNested(depth).sortedExternalEdges(h, delta_depth);
      System.out.println("h: " + h + "; neigs: " + actualRes);
      toAladinDraw(depth + delta_depth, actualRes);
      System.out.println("\n---------");
    }*/
  }
  
  private static final void checkEquals(FlatHashList neig, long[] expected) {
    neig.sortByHashAsc();
    final long[] a = new long[neig.size()];
    neig.arraycopy(0, a, 0, neig.size());
    assert Arrays.equals(a, expected) : "Array: " + Arrays.toString(a);
  }
  
  private void toAladinDraw(final int depth, FlatHashList list) {
    // NOT STANDARD (SINCE WE DO NOT USE INTERVALS!!)
    System.out.format("draw moc %d/", depth);
    FlatHashIterator it = list.iterator();
    for (int i = 0; it.hasNext(); i++) {
      System.out.format("%d,", it.next());
    }
  }
  
  public static final void main(String[] args) {
    final HealpixNestedNeighbourSelectorTest t = new HealpixNestedNeighbourSelectorTest();
    t.externalEdgeTest4();
    t.unsortedExternalEdgeTest4();
    
    
  } 

}
