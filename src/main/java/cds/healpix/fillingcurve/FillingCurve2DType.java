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

package cds.healpix.fillingcurve;

public enum FillingCurve2DType {
  
  Z_ORDER_OR() {
    @Override
    public final FillingCurve2D get(final int depth) {
      if (depth == 0) {
        return ZOrderCurve2DImpls.EMPTY;
      } else if (depth <= 8) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_BYTE;
      } else if (depth <= 16) {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_SHORT;
      } else {
        return ZOrderCurve2DImpls.ZOC_VMSB_OR_INT;
      }
    }
  },
  Z_ORDER_XOR() {
    @Override
    public final FillingCurve2D get(final int depth) {
      if (depth == 0) {
        return ZOrderCurve2DImpls.EMPTY;
      } else if (depth <= 8) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_BYTE;
      } else if (depth <= 16) {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_SHORT;
      } else {
        return ZOrderCurve2DImpls.ZOC_VMSB_XOR_INT;
      }
    }
  },
  Z_ORDER_LUPT() {
    @Override
    public final FillingCurve2D get(final int depth) {
      if (depth == 0) {
        return ZOrderCurve2DImpls.EMPTY;
      } else if (depth <= 8) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_BYTE;
      } else if (depth <= 16) {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_SHORT;
      } else {
        return ZOrderCurve2DImpls.ZOC_VMSB_LOOKUP_INT;
      }
    }
  };
  
  public abstract FillingCurve2D get(final int depth);
  
}
