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

import java.util.Iterator;

import cds.healpix.AbstractHealpixNestedMOC.CurrentCellAccessor;

/**
 * WARNING: IN DEVELOPPEMTENT, NOT YET VISIBLE IN THE JAVADOC
 * TODO: CONTINUE DEV
 * 
 * @author F.-X. Pineau
 *
 */
final class HealpixNestedIntValuedMOC
    extends AbstractHealpixNestedMOC<HealpixNestedIntValuedMOC.CurrentCellWithValAccessor>  {

  
  public static interface CurrentCellWithValAccessor extends CurrentCellAccessor {
    int getAuxValue();
  }
  
  private final class IterVal extends Iter implements CurrentCellWithValAccessor {

    private int auxVal = 0;
    
    @Override
    public final int getAuxValue() {
      return auxVal;
    }
    @Override
    protected final CurrentCellWithValAccessor returnThis() {
      return this;
    }
    @Override
    protected final void calledInNext(final int i) {
      this.auxVal = values[i];
    }
    @Override
    public final String toString() {
      return String.format("d: %d; h: %d; auxVal: %s; raw: %d", getDepth(), getHash(), getAuxValue(), getMOCEncodedHash());
    }
    @Override
    public void remove() {
     throw new UnsupportedOperationException();
    }
  }
  
  
  // protected final int depthMax;
  // protected final long[] cells; // can be int[] or long[]
  protected final int[] values; // can be byte[], short[], int[], long[], float[], double[]
  // protected final int to;

  protected HealpixNestedIntValuedMOC(final int mocDepth, final long[] cells, final int[] values) {
    this(mocDepth, cells, values, cells.length);
  }
  
  protected HealpixNestedIntValuedMOC(final int mocDepth, final long[] cells, final int[] values,
      final int toIndex) {
    super(mocDepth, cells, toIndex);
    this.values = values;
  }
  
  @Override
  public Iterator<CurrentCellWithValAccessor> iterator() {
    return new IterVal();
  }
  
 /* public static interface AuxValueMerger {
    int merger(int auxValueLeft, int auxValueRight);
  }
  public static interface AuxValueFilter {
    int accept(int mergedAuxValue);
  }
  */
  
  /*HealpixNestedIntValuedMOC innerJoin(HealpixNestedIntValuedMOC rightMoc, AuxValueMerger merger, AuxValueFilter filter); // intersection
  HealpixNestedIntValuedMOC leftJoin(HealpixNestedIntValuedMOC rightMoc, AuxValueMerger merger, AuxValueFilter filter);
  HealpixNestedIntValuedMOC leftMinusRight(HealpixNestedIntValuedMOC rightMoc, AuxValueMerger merger, AuxValueFilter filter);
  HealpixNestedIntValuedMOC rightJoin(HealpixNestedIntValuedMOC rightMoc, AuxValueMerger merger, AuxValueFilter filter);
  HealpixNestedIntValuedMOC rightMinusLeft(HealpixNestedIntValuedMOC rightMoc, AuxValueMerger merger, AuxValueFilter filter);
  HealpixNestedIntValuedMOC fullJoin(HealpixNestedIntValuedMOC rightMoc, AuxValueMerger merger, AuxValueFilter filter);  // union
  HealpixNestedIntValuedMOC symDiff(HealpixNestedIntValuedMOC rightMoc, AuxValueMerger merger, AuxValueFilter filter);*/
  
  
  
  // TAKE CARE OF TYPES ONLY WHEN WRITTING/READDING ON DISK
}
