package cds.healpix;

import static org.junit.Assert.assertEquals;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import cds.healpix.NestedEllipticalConeComputerApprox.Mode;

public class NestedEllipticalConeComputerApproxTest {

  
  
  
  @Test
  public void test1() {
    final double aRad = Math.toRadians(14.93);
    final double bRad = Math.toRadians(4.93);
    final double paRad = Math.toRadians(75.0);
    final double lonRad = Math.toRadians(36.80105218);
    final double latRad = Math.toRadians(56.78028536);
    final int depth = 3;
    // draw circle(36.80105218, 56.78028536, 14.93)
    final long[] expectedRes = new long[]{27L, 30L, 39L, 43L, 44L, 45L, 46L, 47L, 48L, 49L, 50L, 
        51L, 52L, 54L, 56L, 57L};
    
    final NestedEllipticalConeComputerApprox cp = new NestedEllipticalConeComputerApprox(
        aRad, bRad,paRad, Healpix.getNested(depth));
    final HealpixNestedBMOC bmoc = cp.overlapping(lonRad, latRad, Mode.OVERLAPPING_CELLS);
    int i = 0;
    //toAladinDraw(depth, bmoc);
    for (final HealpixNestedBMOC.CurrentValueAccessor cell : bmoc) {
      // System.out.println(cell);
      assertEquals(expectedRes[i++], cell.getHash());
    }
    //assertEquals(expectedRes.length, i);
  }
  
  
  private void toAladinDraw (final int depth, HealpixNestedBMOC bmoc) {
    // NOT STANDARD (SINCE WE DO NOT USE INTERVALS!!)
    System.out.format("draw moc %d/", depth);
    FlatHashIterator it = bmoc.flatHashIterator();
    for (int i = 0; it.hasNext(); i++) {
      System.out.format("%d,", it.next());
    }
  }
  
}
