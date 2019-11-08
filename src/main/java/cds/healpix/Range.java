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
}
