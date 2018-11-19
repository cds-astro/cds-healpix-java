package cds.healpix;

/**
 * Interface use for interpolation, i.e. to compute a value (an image pixel for example) from the 
 * values of HEALPix cells.
 * 
 * @author F.-X. Pineau
 *
 */
public interface Interpolator extends HierarchyItem {


  /**
   * From a position on the sky, returns a list of HEALPix cells and associated weights
   * typically depending on the distances between the position and the cell centers.
   * 
   * @param lonRad longitude in radians, must support reasonably large positive and negative values
   *         producing accurate results with a naive range reduction like modulo 2*pi
   *         (i.e. without having to resort on Cody-Waite or Payne Hanek range reduction).
   * @param latRad latitude in [-pi/2, pi/2] radians
   * @param resultHashs array used to store the list of hash to be used for the interpolation
   * @param resultWeights array used to store the weights associated with each hash
   * @return the number of elements to be read in the results arrays
   */
  int interpolate(double lon, double lat,  long[] resultHashs, double[] resultWeights);

  
}
