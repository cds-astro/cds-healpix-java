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

import static cds.healpix.common.math.HackersDelight.SIGN_BIT_MASK_L;
import static cds.healpix.common.math.HackersDelight.fromBits;
import static cds.healpix.common.math.HackersDelight.toBits;

/**
 * Utility class simply defining enums.
 *
 * @author F.-X. Pineau
 */
public final class CompassPoint {

  /*interface ActionBindedToMainWind {
    void northAction();
    void northeastAction();
    void eastAction();
    void southeastAction();
    void southAction();
    void southWestAction();
    void westAction();
    void northwestAction();
  }*/
  
  // Height principal winds (cardinal + ordinal points)
  // ^ SouthToWest axis
  // | _W(0,2) NW(1,2) _N(2,2)
  // | SW(0,1) _C(1,1) NE(2,1)
  // | _S(0,0) SE(1,0) _E(2,0)
  // -----------------------> SouthToEast axis
  /**
   * Enum defining the compass main wind points.
   * 
   * @author F.-X. Pineau
   *
   */
  public enum MainWind {
    /** North. */
    N (2, 2, Cardinal.N) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return northEastValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return northWestValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return northEastValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return northWestValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return northValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return northValue;
      }
    },
    /** Northeast. */
    NE(2, 1, Ordinal.NE) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return northEastValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return centralValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return northEastValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return centralValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return northeastValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return northeastValue;
      }
    },
    /** East. */
    E (2, 0, Cardinal.E) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return northEastValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return southEastValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return northEastValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return southEastValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return eastValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return eastValue;
      }
    },
    /** Southeast. */
    SE(1, 0, Ordinal.SE) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return centralValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return southEastValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return centralValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return southEastValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return southeastValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return southeastValue;
      }
    },
    /** South. */
    S (0, 0, Cardinal.S) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return southWestValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return southEastValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return southWestValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return southEastValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return southValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return southValue;
      }
    },
    /** Southwest. */
    SW(0, 1, Ordinal.SW) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return southWestValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return centralValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return southWestValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return centralValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return southwestValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return southwestValue;
      }
    }, 
    /** West. */
    W (0, 2, Cardinal.W) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return southWestValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return northWestValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return southWestValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return northWestValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return westValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return westValue;
      }
    },
    /** Northwest. */
    NW(1, 2, Ordinal.NW) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return centralValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return northWestValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return centralValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return northWestValue;
      }
      @Override
      
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return northwestValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return northwestValue;
      }
    },
    /** Center. */
    C(1, 1) {
      @Override
      final int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue) {
        return centralValue;
      }
      @Override
      final int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue) {
        return centralValue;
      }
      @Override
      final long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue) {
        return centralValue;
      }
      @Override
      final long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue) {
        return centralValue;
      }
      @Override
      int pickRightIntValue(int northwestValue, int northValue, int northeastValue,
          int westValue, int centralValue, int eastValue,
          int southwestValue, int southValue, int southeastValue) {
        return centralValue;
      }
      @Override
      long pickRightLongValue(long northwestValue, long northValue, long northeastValue,
          long westValue, long centralValue, long eastValue,
          long southwestValue, long southValue, long southeastValue) {
        return centralValue;
      }
    };

    private static final MainWind[] MAIN_WINDS = new MainWind[9];
    static {
      for (final MainWind mainWind : MainWind.values()) {
        MAIN_WINDS[mainWind.index] = mainWind;
      }
      N.setOppositeDirection(S);
      S.setOppositeDirection(N);
      E.setOppositeDirection(W);
      W.setOppositeDirection(E);
      NE.setOppositeDirection(SW);
      SW.setOppositeDirection(NE);
      NW.setOppositeDirection(SE);
      SE.setOppositeDirection(NW);
      C.setOppositeDirection(C);
    }

    private MainWind oppositeDirection;
    private Ordinal ordinal;
    private Cardinal cardinal;
    private final int index;
    private final int iSW;
    private final int iSE;
    private final int offsetSE;
    private final int offsetSW;
    
    private MainWind(int indexSE, int indexSW) {
      this(indexSE, indexSW, null, null);
    }
    private MainWind(int indexSE, int indexSW, Ordinal ordinal) {
      this(indexSE, indexSW, ordinal, null);
    }
    private MainWind(int indexSE, int indexSW, Cardinal cardinal) {
      this(indexSE, indexSW, null, cardinal);
    }
    private MainWind(int indexSE, int indexSW, Ordinal ordinal, Cardinal cardinal) {
      this.iSE = indexSE;
      this.iSW = indexSW;
      this.index = computeIndex(this.iSE, this.iSW);
      this.offsetSE = this.iSE - 1;
      this.offsetSW = this.iSW - 1;
    }
    private void setOppositeDirection(final MainWind oppositeDirection) {
      this.oppositeDirection = oppositeDirection;
    }
    final int getIndex() {
      return this.index;
    }
    /**
     * Returns the opposite direction.
     * @return the opposite direction.
     */
    public MainWind getOppositeDirection() {
      return this.oppositeDirection;
    }
    /**
     * Returns {@code true} is this main wind direction is an ordinal point.
     * @return {@code true} is this main wind direction is an ordinal point.
     */
    public boolean isOrdinal() {
      return this.ordinal != null;
    }
    /**
     * Returns the equivalent ordinal point in the {@link Ordinal} enum
     * (throws an IllegalArgumentException if the direction is not an ordinal point).
     * @return the equivalent ordinal point in the {@link Ordinal} enum.
     */
    public Ordinal toOrdinal() {
      if (!this.isOrdinal()) {
        throw new IllegalArgumentException("Main wind " + this + " is not an ordinal point!");
      }
      return this.ordinal;
    }
    /**
     * Returns {@code true} is this main wind direction is a cardinal point.
     * @return {@code true} is this main wind direction is a cardinal point.
     */
    public boolean isCardinal() {
      return this.cardinal != null;
    }
    /**
     * Returns the equivalent cardinal point in the {@link Cardinal} enum.
     * (throws an IllegalArgumentException if the direction is not an cardinal point.)
     * @return the equivalent cardinal point in the {@link Cardinal} enum.
     */
    public Cardinal toCardinal() {
      if (!this.isCardinal()) {
        throw new IllegalArgumentException("Main wind " + this + " is not a cardinal point!");
      }
      return this.cardinal;
    }
    final int getOffsetSE() {
      return this.offsetSE;
    }
    final int getOffsetSW() {
      return this.offsetSW;
    }
    /**
     * Equivalent to MainWind.values().length.
     * We add this method for performances issues since "value()" returns a newly cloned array.
     * @return the size of the enum, i.e. the number of elements it contains.
     */
    static final int size() {
      return MAIN_WINDS.length;
    }
    /**
     * Returns a compass direction from its position in a 3x3 matrix having the center at the
     * coordinate (1, 1).
     * @param indexSE index along the south to east axis, from 0 (South) to 2 (East).
     * @param indexSW index along the south to west axis, from 0 (South) to 2 (West).
     * @return the compass point a at the given coordinate in a 3x3 matrix.
     */
    static final MainWind getFromCoo(int indexSE, int indexSW) {
      return MAIN_WINDS[computeIndex(indexSE, indexSW)];
    }
    static final MainWind getFromOffset(int offsetSE, int offsetSW) {
      return MAIN_WINDS[computeIndex(++offsetSE, ++offsetSW)];
    }
    static final MainWind getFromIndex(final int i) {
      return MAIN_WINDS[i];
    }
    private static final int computeIndex(int indexSE, int indexSW) {
      // We use (indexSW << 1) + indexSW instead of indexSW * 3 to avoid a multiplication 
      return (indexSW << 1) + indexSW + indexSE;
    }

    abstract int pickRightSouthToEastIntValue(int southWestValue, int centralValue, int northEastValue);
    abstract int pickRightSouthToWestIntValue(int southEastValue, int centralValue, int northWestValue);
    abstract long pickRightSouthToEastLongValue(long southWestValue, long centralValue, long northEastValue);
    abstract long pickRightSouthToWestLongValue(long southEastValue, long centralValue, long northWestValue);
    
    abstract int pickRightIntValue(
        int northwestValue, int   northValue, int northeastValue,
        int      westValue, int centralValue, int      eastValue,
        int southwestValue, int   southValue, int southeastValue);
    abstract long pickRightLongValue(
        long northwestValue, long   northValue, long northeastValue,
        long      westValue, long centralValue, long      eastValue,
        long southwestValue, long   southValue, long southeastValue);

    //abstract <T> T pickRightSouthNorthValue(T southValue, T centralValue, T northValue);
    // abstract void trigger(ActionBindedToMainWind action);
  }

  
  // Ordinal (or intercardinal) points
  /**
   * Enum defining the compass ordinal points.
   *
   * @author F.-X. Pineau
   *
   */
  public enum Ordinal {
    /** Northeast. */
    NE {
      @Override
      void orderedInternalEdge(final NeighbourSelector neigSelect, final long hash,
          final int deltaDepth, final FlatHashList result) {
        neigSelect.sortedInternalEdgeNE(hash, deltaDepth, result);
      }
    },
    /** Southeast. */
    SE {
      @Override
      void orderedInternalEdge(final NeighbourSelector neigSelect, final long hash,
          final int deltaDepth, final FlatHashList result) {
        neigSelect.sortedInternalEdgeSE(hash, deltaDepth, result);
      }
    },
    /** Southwest. */
    SW {
      @Override
      void orderedInternalEdge(final NeighbourSelector neigSelect, final long hash,
          final int deltaDepth, final FlatHashList result) {
        neigSelect.sortedInternalEdgeSW(hash, deltaDepth, result);
      }
    },
    /** Northwest. */
    NW {
      @Override
      void orderedInternalEdge(final NeighbourSelector neigSelect, final long hash,
          final int deltaDepth, final FlatHashList result) {
        neigSelect.sortedInternalEdgeNW(hash, deltaDepth, result);
      }
    };  
    
    abstract void orderedInternalEdge(final NeighbourSelector neigSelect, final long hash,
        final int deltaDepth, final FlatHashList result);

  }
  
  /**
   * Enul defining the compass cardinal points.
   * @author F.-X. Pineau
   *
   */
  public enum Cardinal {
    /** North. */
    N(0,  0,  1) { // North
      @Override
      long internalCorner(final NeighbourSelector neigSelect, final long hash,
          final int toEdgeDeltaDepth) {
        return neigSelect.internalCornerN(hash, toEdgeDeltaDepth);
      }
    },
    /** East. */
    E(1,  1,  0) {
      @Override 
      long internalCorner(final NeighbourSelector neigSelect, final long hash,
          final int toEdgeDeltaDepth) {
        return neigSelect.internalCornerE(hash, toEdgeDeltaDepth);
      }
    },
    /** South. */
    S(2,  0, -1) {
      @Override
      long internalCorner(final NeighbourSelector neigSelect, final long hash,
          final int toEdgeDeltaDepth) {
        return neigSelect.internalCornerS(hash, toEdgeDeltaDepth);
      }
    },
    /** West. */
    W(3, -1,  0) {
      @Override
      long internalCorner(final NeighbourSelector neigSelect, final long hash,
          final int toEdgeDeltaDepth) {
        return neigSelect.internalCornerW(hash, toEdgeDeltaDepth);
      }
    };
    
    private static final Cardinal[] CARDINAL_POINTS = new Cardinal[4];
    static {
      for (final Cardinal cardinalPoint : Cardinal.values()) {
        CARDINAL_POINTS[cardinalPoint.index] = cardinalPoint;
      }
    }
    
    private final int index;
    private final int xOffset;
    private final int yOffset;
    private final long xOffsetSign;
    private final long yOffsetSign;
    private final long xOffsetMak;
    private final long yOffsetMak;
    
    private Cardinal(final int index, final int xAxisOffset, final int yAxisOffset) {
      this.index = index;
      this.xOffset = xAxisOffset;    assert xOffset == -1 || xOffset == 0 || xOffset == 1;
      this.yOffset = yAxisOffset;    assert yOffset == -1 || yOffset == 0 || yOffset == 1;
      // Internal derived variables
      this.xOffsetSign = xAxisOffset >= 0 ? 0 : SIGN_BIT_MASK_L;
      this.xOffsetMak = xAxisOffset == 0 ? 0 : -1;
      this.yOffsetSign = yAxisOffset >= 0 ? 0 : SIGN_BIT_MASK_L;
      this.yOffsetMak = yAxisOffset == 0 ? 0 : -1;
    }
    /**
     * Returns the next cardinal point in the clockwise direction.
     * @return the next cardinal point in the clockwise direction.
     */
    public Cardinal nextClockwise() {
      return CARDINAL_POINTS[(this.index + 1) & 3];
    }
    /**
     * Returns the next cardinal point in the counter-clockwise direction.
     * @return the next cardinal point in the counter-clockwise direction.
     */
    public Cardinal nextCounterClockwise() {
      return CARDINAL_POINTS[(this.index - 1) & 3];
    }
    
    final int getOffsetX() {
      return this.xOffset;
    }
    
    final int getOffsetY() {
      return this.yOffset;
    }
    
    final double timeXOffset(final double absValue) {
      return fromBits((toBits(absValue) | xOffsetSign) & xOffsetMak);
    }

    final double timeYOffset(final double absValue) {
      return fromBits((toBits(absValue) | yOffsetSign) & yOffsetMak);
    }
    abstract long internalCorner(final NeighbourSelector neigSelect, final long hash,
        final int toEdgeDeltaDepth);
  }
  
  /**  Prevents from instantiation. */
  private CompassPoint() { }

}
