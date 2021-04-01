package cds.healpix;

/**
 * Very generic Cell class.
 * The dimension (e.g. HEALPix cell = dim 2; Time cell = dim 1) is not given here.
 * 
 * @author F.-X. Pineau
 *
 */
public class Cell {

    /** Cell depth. */
	public final long depth;
	/** Cell number (hash value). */
	public long hash;
	  
	public Cell(final int depth, final long hash) {
	    this.depth = depth;
	    this.hash = hash;
	}

}
