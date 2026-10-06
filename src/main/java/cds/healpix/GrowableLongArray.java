package cds.healpix;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * This class is used when the number of Hash returned by a function is not necessarily
 * known in advance and may grow.
 * 
 * @author F.-X. Pineau
 *
 */
public final class GrowableLongArray {
    public static final Logger LOGGER = Logger.getLogger( NestedSmallCellApproxedMethod.class.getPackage().getName());
    private long[] array;
    private int cursor;
    public GrowableLongArray(int capacity) {
      assert capacity > 1; // else capacity + (capacity)/2 always returns 1
      this.array = new long[capacity];
      this.cursor = 0;
    }
    public long[] getArray() { return this.array; }
    public int getCursor() { return this.cursor; }
    public final void add(long value) {
      // Mark Taylor suggested to catch the ArrayIndexOutOfBoundsException (I like the idea because 
      // it resorts on the built-in bound check, so we do not have to add an extra test), see 
      //   https://github.com/cds-astro/cds-healpix-java/issues/15
      // I wonder why an explicit test is used in the Java API ArrayList, see e.g.
      //   https://hg.openjdk.java.net/jdk8/jdk8/jdk/file/tip/src/share/classes/java/util/ArrayList.java
      // I guess that it is for better performances when frequent need to make the array grow
      // (but in our case the operation is supposed to be infrequent).
      // Is the compiler/jit smart enough not to make the test (test + bound check) twice?
      // I so far leave as it is: ideally I should have checked performances with both solutions
      // (difference probably negligible, but to be checked!).
      if (this.cursor == this.array.length) {
        // There is no debug() method (no DEBUG Level) in the java.util.logging.Logger.
        // One can use fine/finer/finest instead.
        // Java default output is INFO, see https://docs.oracle.com/cd/E17277_02/html/GettingStartedGuide/managelogging.html) 
        // So this will not show up by default (but can be activated if needed)
        LOGGER.finer("Had to grow unpacked moc size!");
        // New size = old size + old size / 2 (same code as in Java API ArrayList)
        this.array = Arrays.copyOf(this.array, this.array.length + (this.array.length >> 1));
      }
      this.array[this.cursor++] = value;
    }
  }