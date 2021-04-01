package cds.healpix;

import java.util.Iterator;

/**
 * Range iterator to cell iterator decorator.
 * We could have made this decorator more generic by replacing both the 
 * Range and the Cell classes by an interface, adding a object building a Cell from a Range
 * plus the cell depth and hash.
 * 
 * 
 * @author F.-X Pineau
 */
public final class ToCellItDecorator implements Iterator<Cell> {

  private final Iterator<Range> it;
  private int shift;
  private int absoluteDephtMax;
  
  private long l = 0, h = 0;
  private Cell current;
  private Cell next;
	
  /**
   * 
   * @param dim dimension of the z-order curve. For HEALPix, dim = 2, for Time, dim = 1
   * @param absoluteDephtMax for HEALPix, absoluteDephtMax = 29
   * @param decorated
   */
  public ToCellItDecorator(final int dim, final int absoluteDephtMax, final Iterator<Range> decorated) {
	this.it = decorated;
	this.shift = (dim - 1);
	this.absoluteDephtMax = absoluteDephtMax;
    if (this.it.hasNext()) {
      final Range r = this.it.next();
      this.l = r.from;
      this.h = r.to;
      this.next = uncheckedNext();
	}
	this.nextCell();
  }

  private final Cell uncheckedNext() {
	long len = this.h - this.l;
    int ddMaxFromLen = (63 - Long.numberOfLeadingZeros(len)) >> this.shift;
    int ddMaxFromLow = Long.numberOfTrailingZeros(l) >> this.shift;
    int dd = Math.min(this.absoluteDephtMax, Math.min(ddMaxFromLen, ddMaxFromLow));
    int ddShift = dd << this.shift;
    final Cell res = new Cell(this.absoluteDephtMax - dd, l >> ddShift);
    this.l += 1L << ddShift;
    return res;
  }
  
  private final void nextCell() {
    this.current = this.next;
	if (this.l < this.h) {
      this.next = uncheckedNext();
    } else if (this.hasNext()) {
      final Range r = this.it.next();
      this.l = r.from;
      this.h = r.to;
      this.next = uncheckedNext();
    } else {
      this.next = null;	
    }
  }

  @Override
  public boolean hasNext() {
	return this.current != null;
  }

  @Override
  public Cell next() {
	final Cell r = this.current;
    this.nextCell();
    return r;
  }

}
