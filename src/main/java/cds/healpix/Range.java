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
   *  
   * @param sink
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
      l += 1 << twiceDd;
    } while (l < h);
  }
  
}
