package cds.healpix;

/**
 * level 0:
 * - 3 long storing one bit by cell of level 2
 * level 1:
 * - if orderMax = 3 => 1 long storing 16 sub-cells of level 3 (16 * 4     = 64) 
 * - if orderMax = 4 => 1 long storing  4 sub-cells of level 4 ( 4 * 4*4   = 64)
 * - if orderMax = 5 => 1 long storing  1 sub-cell  of level 5 ( 1 * 4*4*4 = 64)
 * - if orderMax > 5 => 
 *   -  1 long  storing 1 sub-cell of level 4 ( 1 * 4*4*4 = 64)
 *   - +1 short or byte storing the cumulative number of +1 
 * level 2:
 * - if orderMax = 6 => 1 long storing 16 sub-cells of level 6 (16 * 4     = 64) 
 * - if orderMax = 7 => 1 long storing  4 sub-cells of level 7 ( 4 * 4*4   = 64)
 * - if orderMax = 8 => 1 long storing  1 sub-cell  of level 8 ( 1 * 4*4*4 = 64)
 * - if orderMax > 8 => 
 *   -  1 long  storing 1 sub-cell of level 4 ( 1 * 4*4*4 = 64)
 *   - +1 int or short or byte storing the cumulative number of +1
 * level 3:
 *  ...
 * 
 * @author F.-X. Pineau
 *
 */
class HealpixNestedCompactStructuredMOC {
  
  private HealpixNestedCompactStructuredMOC() {
    
  }
  
  private static HealpixNestedCompactStructuredMOC from(final HealpixNestedBMOC bmoc) {
    final int depthMax = bmoc.getDepthMax();
    final int treeDepthMax = (int) ((depthMax) / 3); // = 0 (0-2), 1 (3-5), 2 (6-8), 3 (9-11), 4 (12-14), ...
    
    if (true) {
      
    } else {
      
    }
    
    return null;
  }
  
}
