package cds.healpix;

import java.util.ArrayList;
import java.util.Arrays;
// import java.util.Base64;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;

import cds.healpix.HealpixNestedBMOC.CurrentValueAccessor;

/**
 * 
 * The idea is to use an implicit datastructure, encoding the natural tree traversal on bits with:
 * - 0 bit  => cell not in the MOC, go to the next sibling cell or go up if all siblings already explored
 * - 1 bit  => cell contains information, go deeper
 * - 0000 bits => parent cell is in the MOC, no need to go deeper (ony for cells of depth < depthMax)
 * 
 * @author F.-X. Pineau
 *
 */
public class HealpixNestedUltraCompactMOC {
  
  private static final class CustomBitSet {
    
    private BitSet bs;
    private int capacity;
    private int i = 0;
    
    private CustomBitSet(final int mocSize) {
      this.bs = new BitSet(64 * mocSize);
      this.capacity = this.bs.size();
    }
    
    private void add(boolean bitValue) {
      if (this.i >= this.capacity) {
        final int newCapacity = this.capacity + (int) (0.5 * this.capacity);
// System.out.println("new capacity: " + newCapacity);
        final BitSet newBitSet = new BitSet(newCapacity);
        newBitSet.or(this.bs);
        this.bs = newBitSet;
        this.capacity = newCapacity;
      }
      this.bs.set(this.i++, bitValue);
    }
    
    private byte[] toByteArray() {
      // Replace to stay compatible with Java6 
      // return Arrays.copyOf(this.bs.toByteArray(), (int) ((i + 7) / 8));
      byte[] bytes = new byte[(int) ((i + 7) / 8)];
      for (int i = 0; i < this.i; i++) {
        if (this.bs.get(i)) {
          bytes[bytes.length - i / 8 - 1] |= 1 << (i & 7);
        }
      }
      return bytes;
    }
  }
  
  
  public static byte[] compress(final HealpixNestedBMOC moc) {
    final CustomBitSet bits = new CustomBitSet(moc.size());
    final int depthMax = moc.getDepthMax();
    final Iterator<CurrentValueAccessor> it = moc.iterator();
    CurrentValueAccessor curr;
    int depth = 0, currDepth = 0;
    long hash = 0, currHash = 0;
    // start, go down
    if (it.hasNext()) {
      curr = it.next();
      currDepth = curr.getDepth();
      currHash = curr.getHash();
      // go down to currHash
      for (int dd = currDepth - depth; dd >= 0; hash <<= 2, dd--) {
        for (long targetHash = currHash >>> (dd << 1); hash < targetHash; hash++) {
          bits.add(false);
        }
        bits.add(true);
      }
      if (currDepth != depthMax) {
        bits.add(false);
        bits.add(false);
        bits.add(false);
        bits.add(false);
      }
      hash = currHash;
      depth = currDepth;
    }
    // middle, go up and down
    for (curr = it.next(), currDepth = curr.getDepth(), currHash = curr.getHash()
        ; it.hasNext()
        ; hash = currHash, depth = currDepth,
         curr = it.next(), currDepth = curr.getDepth(), currHash = curr.getHash()) {
      // go up (if needed)!
      int dd = currDepth - depth;
      long currHashAtPrevDepth;
      if (dd < 0) { // case previous hash deeper that current hash
        int nBits = (-dd) << 1;
        currHashAtPrevDepth = currHash << nBits;
      } else {      // case current hash deeper that previous hash, need to go up?
        currHashAtPrevDepth = currHash >>> (dd << 1);
      }
      dd = ((63 - Long.numberOfLeadingZeros(hash ^ currHashAtPrevDepth)) >> 1);
      if (dd > depth) {
        dd = depth;
      }
      depth -= dd;
      // - go up to depth common depth
      for (; dd > 0; hash >>= 2, dd--) {
        for (; (hash & 3L) < 3L; hash++) {
          bits.add(false);
        }
      }
      hash++;
      // go down!
      for (dd = currDepth - depth; dd >= 0; hash <<= 2, dd--) {
        for (long targetHash = currHash >>> (dd << 1); hash < targetHash; hash++) {
          bits.add(false);
        }
        bits.add(true);
      }
      if (currDepth != depthMax) {
        bits.add(false);
        bits.add(false);
        bits.add(false);
        bits.add(false);
      }
    }
    // end, go up
    // - go up to depth 0
    for (; depth > 0; hash >>= 2, depth--) {
      for (int k = (((int) hash) & 3); k < 3; k++) {
        bits.add(false);
      }
    }
    hash++;
    // - complete till base cell 11
    for (; hash < 12; hash++) {
      bits.add(false);
    }
    return bits.toByteArray();
  }
 
  private static final BitSet bitSetFromBytes(byte[] compressedMoc) {
    final BitSet bs = new BitSet(8 * compressedMoc.length);
    for (int i = 0; i < bs.size(); i++) {
      bs.set(i,((compressedMoc[i / 8] & (1 << (i & 7))) != 0));
    }
    return bs;
  }
  
  public static HealpixNestedBMOC decompress(int depthMax, byte[] compressedMoc) {
    final BitSet bs = bitSetFromBytes(compressedMoc); // BitSet.valueOf(compressedMoc); CHANGED TO BE COMPATIBLE WITH JAVA6
    final List<Long> res = new ArrayList<Long>(compressedMoc.length);
    int i = 0;
    int depth = 0;
    long hash = 0;
    do {
      boolean bit = bs.get(i);
      ++i;
      if (bit) { // = 1
        if (depth == depthMax) {
          res.add(HealpixNestedBMOC.buildValue(depth, hash, true, depthMax));
          for (; (hash & 3L) == 3L && depth > 0; hash >>>= 2, depth--); // go up if needed
          ++hash;
        } else {
          // go down of 1 level
          hash <<= 2;
          depth++;
        }
      } else { // = 0
        if (depth == 0) {
          if (hash == 12) { break; }
          ++hash;
        } else {
          // Case 0000 => hash of deph = depth - 1 fully in MOC
          if ((hash & 3L) == 0L && !bs.get(i) && !bs.get(i + 1) && !bs.get(i + 2)) {
            --depth;
            hash >>>= 2;
            i += 3;
            res.add(HealpixNestedBMOC.buildValue(depth, hash, true, depthMax));
          } 
          // Case hash = xxx...11
          for (; (hash & 3L) == 3L && depth > 0; hash >>>= 2, depth--); // go up if needed
          ++hash;
        }
      }
    } while (true);
    final int size = res.size();
    final long[] mocCells = new long[size];
    for (int k = 0; k < size; k++) {
      mocCells[k] = res.get(k);
    }
    return HealpixNestedBMOC.createUnsafe(depthMax, mocCells);
  }
  
  
  /* REMOVED To STAY COMPATIBLE WITH JAVA6
  public static String compressB64(byte[] bytes) {
    return Base64.getEncoder().encodeToString(bytes);
  }
  
  public static byte[] decompressB64(final String b64Encoded) {
    return Base64.getDecoder().decode(b64Encoded);
  }*/
  
}
