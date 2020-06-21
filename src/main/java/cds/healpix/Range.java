package cds.healpix;

/**
 * Defines an Hash Range at the deeper depth (i.e. 29).
 * @author F.-X. Pineau
 *
 */
public class Range {
  
  /** Range lower bound (inclusive). */
  public final long from;
  /** Range lower bound (exclusive). */
  public long to;
  
  public Range(final long from, final long to) {
    this.from = from;
    this.to = to;
  }
  
  /**
   * Transforms this range in a list of cells that are added to the given {@code sink}.
   * IMPORTANT: the order in which the cells are added follows the natural Z-order curve order!
   * @param sink object receiving each cell
   */
  public void toCells(final CellSink sink) {
    long l = this.from;
    long h = this.to;
    do {
      long len = h - l;
      assert len > 0;
      int ddMaxFromLen = (63 - Long.numberOfLeadingZeros(len)) >> 1;
      int ddMaxFromLow = Long.numberOfTrailingZeros(l) >> 1;
      int dd = Math.min(29, Math.min(ddMaxFromLen, ddMaxFromLow));
      int twiceDd = dd << 1;
      sink.push(29 - dd, l >> twiceDd);
      l += 1L << twiceDd;
    } while (l < h);
  }
  
  /**
   * Same as {@code toCells} but with additional informations which are (see parameter list).
   * This version may have better performances since a large part of the cells in a MOC are
   * at the deepest MOC order.
   * @param sink object receiving each cell
   * @param depthMax the depth of the lower possible cell order in the Range (i.e. the MOC order)
   * @param twiceDD {@code (29 - depthMax) << 1}, provided not ot have to recompute it
   * @param rangeLenMin {@code 1L << twiceDD}, provided not to have to recompute it
   * @param mask {@code 3L << twiceDD}, provided not to have to recompute it
   */
  public void toCellsWithKnowledge(final CellSink sink,
      int depthMax, int twiceDD, long rangeLenMin, long mask) {
    long l = this.from;
    long h = this.to;
    do {
      long len = h - l;
      assert len > 0;
      if (len == rangeLenMin || (l & mask) != 0L) {
        // A range of 1 cell at depthMax
        sink.push(depthMax, l >> twiceDD);
        l += rangeLenMin;
      } else {
        int ddMaxFromLen = (63 - Long.numberOfLeadingZeros(len)) >> 1;
        int ddMaxFromLow = Long.numberOfTrailingZeros(l) >> 1;
        int dd = Math.min(29, Math.min(ddMaxFromLen, ddMaxFromLow));
        int twiceDd = dd << 1;
        sink.push(29 - dd, l >> twiceDd);
        l += 1L << twiceDd; 
      }
    } while (l < h);
  }
}
