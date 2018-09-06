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

import static cds.healpix.Healpix.getBestStartingDepth;

/**
 * The idea of this class is to avoid making multiple time the same operations (like selecting
 * the optimal starting depth) in case of fixed radius cross-match.
 * 
 * This is especially made for cross-matches with the folowing ideas:
 *  - The list of cells returned is small (maximum 9 elements)
 *  - we assume that cells are not fully overlaped by the xmatch cone
 *  - we accept false positives (i.e. cell wich are close but do not overlap the xmatch cone
 * 
 * The idea is to have a function which is a quick as possible (no much longer than looking
 * into a HashMap).
 * 
 * @author F.-X. Pineau
 *
 */
public final class HealpixNestedFixedRadiusCone4XMatch {

  private final int hashDepth;
  private final HealpixNested hn;
  private final HashComputer hc;
  private final NeighbourSelector neiSelect;
  private final int twiceDeltaOrder;
  private final double coneRadiusRad;
  private final NeighbourList neigList;
  
  protected HealpixNestedFixedRadiusCone4XMatch(final int hashDepth, final double coneRadiusRad) {
    this.hashDepth = hashDepth;
    final int optimalStartingDepth = getBestStartingDepth(coneRadiusRad);
    if (optimalStartingDepth < this.hashDepth) {
      throw new IllegalArgumentException("Choosen radius too large for the wanted depth!");
    }
    this.hn = Healpix.getNested(optimalStartingDepth);
    this.hc = this.hn.newHashComputer();
    this.neiSelect = this.hn.newNeighbourSelector();
    this.twiceDeltaOrder = (optimalStartingDepth - hashDepth) << 1;
    this.coneRadiusRad = coneRadiusRad;
    this.neigList = new NeighbourList(this.hashDepth);
  }
  
  private HealpixNestedFixedRadiusCone4XMatch(final HealpixNestedFixedRadiusCone4XMatch o) {
    this.hashDepth = o.hashDepth;
    this.hn = o.hn;
    this.hc = this.hn.newHashComputer();
    this.neiSelect = this.hn.newNeighbourSelector();
    this.twiceDeltaOrder = o.twiceDeltaOrder;
    this.coneRadiusRad = o.coneRadiusRad;
    this.neigList = new NeighbourList(this.hashDepth);
  }
  
  /**
   * Returns the radius of the cones, in radians.
   * @return he radius of the cones, in radians.
   */
  public double getRadius() {
    return this.coneRadiusRad;
  }
  
  public int getResultDepth() {
    return this.hashDepth;
  }
  
  /**
   * Fill the given array with the cells overallping the cone of given center
   * (may include false positive). We assume that the size of a cell at the depth of the returned
   * hashs is large compared to the cone radius.
   * @param coneCenterLonRad longitude of the center of the cone, in radians
   * @param coneCenterLatRad latitude of the center of the cone, in radians
   * @param result an array of size al least 9 that will store overlapping cell from index 0,
   *        the content of the array is NOT sorted.
   * @return the number of elements to be read in the list
   */
  public int overlappingCells(final double coneCenterLonRad, final double coneCenterLatRad, long[] result) {
    final long centerHash = this.hc.hash(coneCenterLonRad, coneCenterLatRad);
    this.neiSelect.neighbours(centerHash, this.neigList);
    final long ch = centerHash >>> this.twiceDeltaOrder;
    result[0] = ch;
    int ir = 1;
    for (int i = 0; i < this.neigList.size(); i++) {
      final long h = this.neigList.get(i) >>> this.twiceDeltaOrder;
      if (h != ch && isNotIn(h, result, ir)) {
        result[ir++] = h;
      }
    }
    return ir;
  }
  private static final boolean isNotIn(final long h, final long[] ah, final long alength) {
    for (int i = 1; i < alength; i++) {
      if (h == ah[i]) {
        return false;
      }
    }
    return true;
  }
  
  /**
   * To obtain new instances in case we want to use multi-threading,
   * since an object is possibly not thread-safe.
   * If the objectis thread-safe, the method can simply return {@code this}.
   * @return a new instance of this class.
   */
  public HealpixNestedFixedRadiusCone4XMatch newComputer() {
    return new HealpixNestedFixedRadiusCone4XMatch(this);
  }

  
}
