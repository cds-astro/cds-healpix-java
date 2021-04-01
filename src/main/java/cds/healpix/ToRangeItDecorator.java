package cds.healpix;

import java.util.Iterator;

import cds.healpix.HealpixNestedBMOC.CurrentValueAccessor;

/**
 * WARNING: this is so far dim = 2 specific!
 * @author F.-X Pineau
 *
 */
public class ToRangeItDecorator implements Iterator<Range> {
  
  private final Iterator<CurrentValueAccessor> it;
  private Range current;
  private Range next;


  public ToRangeItDecorator(final Iterator<CurrentValueAccessor> decorated) {
    this.it = decorated;
    if (this.it.hasNext()) {
      this.next = toRange(this.it.next());
    }
    this.nextRange();
  }
  
  private void nextRange() {
    Range r = null;
    while (this.it.hasNext() && tryMerge(this.next, r = toRange(this.it.next()))) {
      r = null; // to avoid case in which last elem merge (=> is set in the this.next variable)
    }
    this.current = this.next;
    this.next = r;
  }
  
  private Range toRange(final CurrentValueAccessor cva) {
    final int depth = cva.getDepth();
    final long hash = cva.getHash();
    final int ddTwice = (29 - depth) << 1;
    return new Range(hash << ddTwice, (hash + 1L) << ddTwice);
  }
  
  /**
   * 
   * @param left
   * @param right
   * @return {@code true} if the right range has been merge with the left range.
   * The result is stored in the left range.
   */
  private boolean tryMerge(Range left, Range right) {
    if (left.to == right.from) {
      left.to = right.to;
      return true;
    } else {
      return false;
    }
  }
  
  @Override
  public boolean hasNext() {
    return this.current != null;
  }

  @Override
  public Range next() {
    final Range r = this.current;
    this.nextRange();
    return r;
  }
  
}
