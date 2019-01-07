// Copyright 2017-2018 - Universite de Strasbourg/CNRS
// The CDS HEALPix library is developped by the Centre de Donnees
// astronomiques de Strasbourgs (CDS) from the following external papers:
//  - [Gorsky2005]     - "HEALPix: A Framework for High-Resolution Discretization and
//                       Fast Analysis of Data Distributed on the Sphere"
//                       http://adsabs.harvard.edu/abs/2005ApJ...622..759G
//  - [Calabretta2004] - "Mapping on the HEALPix grid"
//                       http://adsabs.harvard.edu/abs/2004astro.ph.12607C
//  - [Calabretta2007] - "Mapping on the HEALPix grid"
//                       http://adsabs.harvard.edu/abs/2007MNRAS.381..865C
//  - [Reinecke2015]   - "Efficient data structures for masks on 2D grids"
//                       http://adsabs.harvard.edu/abs/2015A&A...580A.132R
// It is distributed under the terms of the BSD License 2.0
//
// This file is part of the CDS HEALPix library.
//

package cds.healpix;

import java.util.Arrays;
import java.util.Iterator;

import cds.healpix.AbstractHealpixNestedMOC.CurrentCellAccessor;
import cds.healpix.HealpixNestedBMOC.CurrentValueAccessor;

import static cds.healpix.AbstractHealpixNestedMOC.encodeHash4MOC;

/**
 * WARNING: IN DEVELOPPEMTENT, NOT YET VISIBLE IN THE JAVADOC
 * TODO: CONTINUE DEV
 * 
 * For 'AND', 'OR', 'NOT', 'XOR' and 'MINUS' operations, it is probably much faster to work on
 * ranges. This implementation is made at least to test the performances difference.
 * 
 * @author F.-X. Pineau
 *
 */
final class HealpixNestedMOC extends AbstractHealpixNestedMOC<CurrentCellAccessor> {

  private HealpixNestedMOC(final int mocDepth, final long[] mocCells) {
    this(mocDepth, mocCells, mocCells.length);
  }

  private HealpixNestedMOC(final int mocDepth, final long[] mocCells, final int toIndex) {
    super(mocDepth, mocCells, toIndex);
  }

  public HealpixNestedMOC toDeeperDepth(final int newDepth) {
    if (newDepth <= this.depthMax) {
      throw new IllegalArgumentException("The given depth must be higher than the depth max of the MOC");
    }
    final int twiceDepthDiff = (newDepth - this.depthMax) << 1; 
    final long[] mocNew = new long[this.to];
    for (int i = 0; i < this.to; i++) {
      long c = this.cells[i];
      mocNew[i] = c << twiceDepthDiff;
    }
    return new HealpixNestedMOC(newDepth, mocNew);
  }
  
  public HealpixNestedMOC toLowerDepth(final int newDepth) {
    if (newDepth >= this.depthMax) {
      throw new IllegalArgumentException("The given depth lust be lower than the depth max of the MOC");
    }
    final long[] mocNew = new long[this.to];
    int iNew = 0;
    long prevHashAtNewDepth = -1, currHashAtNewDepth;
    for (int i = 0; i < this.to; i++) {
      final long raw = this.cells[i];
      final int depth = getDepth(raw, this.depthMax);
      if (depth <= newDepth) {
        mocNew[iNew++] = raw;
      } else if (prevHashAtNewDepth != 
          (currHashAtNewDepth = getHashFromDepth(raw, depth) >> (depth - newDepth))) {
        mocNew[iNew++] = (prevHashAtNewDepth << 1) | 1L;
        prevHashAtNewDepth = currHashAtNewDepth;
      }
    }
    return newFrom(newDepth, mocNew, iNew);
  }

  /**
   * The depth of the returned MOC is the larger depth between this MOC and the given MOC.
   * @param rightMoc
   * @return
   */
  public HealpixNestedMOC and(final HealpixNestedMOC rightMoc) { // = intersect = inner join
    final int depthMax = Math.max(this.depthMax, rightMoc.depthMax);
    final long[] result = new long[Math.min(this.size(), rightMoc.size())];
    int iResult = 0;
    final Iterator<CurrentCellAccessor> leftIt = this.iterator(), rightIt = rightMoc.iterator();
    // Test for the easy case: one of the two MOC is empty
    if (!leftIt.hasNext() || !rightIt.hasNext()) {
      return new HealpixNestedMOC(depthMax, result, iResult);
    }
    // Too bad: we have to do the job
    CurrentCellAccessor currLeft = leftIt.next(), currRight = rightIt.next();
    int dl =  currLeft.getDepth(), dr = currRight.getDepth();
    long hl =  currLeft.getHash(), hr = currRight.getHash();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    // The code contains ugly repetitions :o/ (I am looking forward to write this in Rust!) 
    for (;;) {
      if (dl < dr) {
        final int dd = dr - dl;
        final long hrAtdl = hr >>> (dd << 1);
        if (hl < hrAtdl) {
          if (!leftIt.hasNext()) { break; }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hl > hrAtdl) {
          if (!rightIt.hasNext()) { break; }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hl == hrAtdl;
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) { break; }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        }
      } else if (dl > dr) {
        final int dd = dl - dr;
        final long hlAtdr = hl >>> (dd << 1);
        if (hlAtdr < hr) {
          if (!leftIt.hasNext()) { break; }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hlAtdr > hr) {
          if (!rightIt.hasNext()) { break; }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hlAtdr == hr;
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) { break; }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        }
      } else {
        assert dl == dr;
        if (hl < hr) {
          if (!leftIt.hasNext()) { break; }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hl > hr) {
          if (!rightIt.hasNext()) { break; }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hl == hr;
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext() || !rightIt.hasNext()) { break; }
          currLeft  =  leftIt.next();  dl =  currLeft.getDepth(); hl =  currLeft.getHash();
          currRight = rightIt.next();  dr = currRight.getDepth(); hr = currRight.getHash();
        }
      }
    }
    return newFrom(depthMax, result, iResult);
  }
  
  public HealpixNestedMOC or(final HealpixNestedMOC rightMoc) { // = union = fullJoin
    final int depthMax = Math.max(this.depthMax, rightMoc.depthMax);
    final Iterator<CurrentCellAccessor> leftIt = this.iterator(), rightIt = rightMoc.iterator();
    // Easy cases: one of the MOC is empty
    if (!leftIt.hasNext()) { // This MOC is empty
      return new HealpixNestedMOC(depthMax, Arrays.copyOf(rightMoc.cells, rightMoc.to));
    } else if (!rightIt.hasNext()) { // The external MOC is empty
      return new HealpixNestedMOC(depthMax, Arrays.copyOf(this.cells, this.to));
    }
    // Too bad, we actually have to do the job
    final long[] result = new long[this.size() + rightMoc.size()]; // Upper bound on the MOC size
    int iResult = 0;
    CurrentCellAccessor currLeft = leftIt.next(), currRight = rightIt.next();
    int dl =  currLeft.getDepth(), dr = currRight.getDepth();
    long hl =  currLeft.getHash(), hr = currRight.getHash();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL
    // The code contains ugly repetitions :o/ (I am looking forward to write this in Rust!) 
    for (;;) {
      if (dl < dr) {
        final int dd = dr - dl;
        final long hrAtdl = hr >>> (dd << 1);
        if (hl < hrAtdl) {
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hrAtdl < hl) {
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hl == hrAtdl;
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          do {
            if (!rightIt.hasNext()) {
              iResult = completeWith(leftIt, result, iResult);
              break;
            }
            currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
          } while (contains(dl, hl, dr, hr));
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        }
      } else if (dl > dr) {
        final int dd = dl - dr;
        final long hlAtdr = hl >>> (dd << 1);
        if (hlAtdr < hr) {
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hr < hlAtdr) {
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hlAtdr == hr;
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          do {
            if (!leftIt.hasNext()) {
              iResult = completeWith(rightIt, result, iResult);
              break;
            }
            currLeft = leftIt.next(); dl = currLeft.getDepth(); hl = currLeft.getHash();
          } while (contains(dr, hr, dl, hl));
          currRight = rightIt.next();  dr = currRight.getDepth();  hr = currRight.getHash();
        }
      } else {
        assert dl == dr;
        if (hl < hr) {
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hr < hl) {
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hl == hr;
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          } else if (!rightIt.hasNext()) { 
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currLeft  =  leftIt.next();  dl =  currLeft.getDepth(); hl =  currLeft.getHash();
          currRight = rightIt.next();  dr = currRight.getDepth(); hr = currRight.getHash();
        }
      }
    }
    return newFrom(depthMax, result, iResult);
  }
  
  private static int completeWith(final Iterator<CurrentCellAccessor> it, final long[] result, int iResult) {
    while (it.hasNext()) {
      result[iResult++] = it.next().getMOCEncodedHash();
    }
    return iResult;
  }
  
  private static final boolean contains(
      final int largeCellDepth, final long largeCellHash,
      final int smallCellDepth, final long smallCellHash) {
    return largeCellDepth <= smallCellDepth 
        && largeCellHash == smallCellHash >>> ((smallCellDepth - largeCellDepth) << 1);
  }
  
  /*private static final boolean contains(final int largeCellDepth, final long largeCellHash,
      final CurrentCellAccessor smallCellAccessor, int smallCellDepth, long smallCellHash) {
    assert smallCellAccessor.getDepth() == smallCellDepth; // ugly I know :o/
    assert smallCellAccessor.getHash() == smallCellHash;   // ugly I know :o/
    return smallCellDepth >= largeCellDepth 
        && largeCellHash == smallCellHash >>> ((smallCellDepth - largeCellDepth) << 1);
  }*/
  
  
  public HealpixNestedMOC xor(final HealpixNestedMOC rightMoc) { // difference = (A or B) - (A and B)
    // if == do bot put
    // else check for overlap
    // - if overlap, decompose (go down) the bigger cell (and keep trace to possibly go up after!) 
    
    final int depthMax = Math.max(this.depthMax, rightMoc.depthMax);
    final Iterator<CurrentCellAccessor> leftIt = this.iterator(), rightIt = rightMoc.iterator();
    // Easy cases: one of the MOC is empty
    if (!leftIt.hasNext()) { // This MOC is empty
      return new HealpixNestedMOC(depthMax, Arrays.copyOf(rightMoc.cells, rightMoc.to));
    } else if (!rightIt.hasNext()) { // The external MOC is empty
      return new HealpixNestedMOC(depthMax, Arrays.copyOf(this.cells, this.to));
    }
    // Too bad, we actually have to do the job
    final long[] result = new long[this.size() + rightMoc.size()]; // Upper bound on the MOC size
    int iResult = 0;
    CurrentCellAccessor currLeft = leftIt.next(), currRight = rightIt.next();
    int dl =  currLeft.getDepth(), dr = currRight.getDepth();
    long hl =  currLeft.getHash(), hr = currRight.getHash();
    // We have 9 cases to take into account:
    // -  3: dL == dR, dL < dR and dR < dL
    // - x3: hL == hR, hL < hR and hR < hL            dr = currRight.getDepth(); hr = currRight.getHash();

    // The code contains ugly repetitions :o/ (I am looking forward to write this in Rust!) 
    for (;;) {
      if (dl < dr) {
        final int dd = dr - dl;
        final long hrAtdl = hr >>> (dd << 1);
        if (hl < hrAtdl) {
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hrAtdl < hl) {
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hl == hrAtdl;
          // go down from hl (at dl) to hr at dr
          long prevh = hl << 2;
          for (int currd = dl + 1, ddm = dd - 1; ddm >= 0; prevh <<= 2, currd++) {
            for (long targetHash = hr >>> (ddm << 1); prevh < targetHash; prevh++) {
              result[iResult++] = encodeHash4MOC(currd, prevh, this.depthMax);
            }
          }
          // for OR, not for XOR: result[iResult++] = currLeft.getMOCEncodedHash();
          for (;;) {
            int prevd = dr;
            prevh = hr;
            if (!rightIt.hasNext()) {
              // - go up to common depth
              for (; dl < prevd && (prevh & 3L) == 3L; prevh >>= 2, prevd--); // go up to first level to be completed
              for (; dl < prevd; prevh >>= 2, prevd--) { // complete going up
                for (; (++prevh & 3L) <= 3L; ) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
              // complete with left
              iResult = completeWith(leftIt, result, iResult);
              break;
            }
            currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
            if (contains(dl, hl, dr, hr)) {
              // go up to common depth
              int ddm = dr - prevd;
              long hAtPrevDepth;
              if (ddm < 0) { // case previous hash deeper that current hash
                int nBits = (-ddm) << 1;
                hAtPrevDepth = hr << nBits;
              } else {      // case current hash deeper that previous hash, need to go up?
                hAtPrevDepth = hr >>> (ddm << 1);
              }
              ddm = ((63 - Long.numberOfLeadingZeros(prevh ^ hAtPrevDepth)) >> 1);
              // - handle case diff is in base cell
              if (ddm > prevd) {
                ddm = prevd;
              }
              // - go up to common depth
              for (; ddm > 0; prevh >>= 2, ddm--, prevd--) {
                for (; (prevh & 3L) < 3L; prevh++) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
              prevh++;
              // go down!
              for (ddm = dr - prevd; ddm >= 0; prevh <<= 2, ddm--) {
                for (long targetHash = hr >>> (ddm << 1); prevh < targetHash; prevh++) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
            } else {
              // - go up to common depth
              for (; dl < prevd && (prevh & 3L) == 3L; prevh >>= 2, prevd--); // go up to first level to be completed
              for (; dl < prevd; prevh >>= 2, prevd--) { // complete going up
                for (; (++prevh & 3L) <= 3L; ) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
            }
          }          
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        }
      } else if (dl > dr) {
        final int dd = dl - dr;
        final long hlAtdr = hl >>> (dd << 1);
        if (hlAtdr < hr) {
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hlAtdr > hr) {
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hlAtdr == hr;
         // go down from hr (at dr) to hl at dl
          long prevh = hr << 2;
          for (int currd = dr + 1, ddm = dd - 1; ddm >= 0; prevh <<= 2, currd++) {
            for (long targetHash = hl >>> (ddm << 1); prevh < targetHash; prevh++) {
              result[iResult++] = encodeHash4MOC(currd, prevh, this.depthMax);
            }
          }
          // for OR, not for XOR: result[iResult++] = currLeft.getMOCEncodedHash();
          for (;;) {
            int prevd = dl;
            prevh = hl;
            if (!leftIt.hasNext()) {
              // - go up to common depth
              for (; dr < prevd && (prevh & 3L) == 3L; prevh >>= 2, prevd--); // go up to first level to be completed
              for (; dr < prevd; prevh >>= 2, prevd--) { // complete going up
                for (; (++prevh & 3L) <= 3L; ) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
              // complete with left
              iResult = completeWith(rightIt, result, iResult);
              break;
            }
            currLeft = leftIt.next(); dr = currLeft.getDepth(); hr = currLeft.getHash();
            if (contains(dr, hr, dl, hl)) {
              // go up to common depth
              int ddm = dl - prevd;
              long hAtPrevDepth;
              if (ddm < 0) { // case previous hash deeper that current hash
                int nBits = (-ddm) << 1;
                hAtPrevDepth = hl << nBits;
              } else {      // case current hash deeper that previous hash, need to go up?
                hAtPrevDepth = hl >>> (ddm << 1);
              }
              ddm = ((63 - Long.numberOfLeadingZeros(prevh ^ hAtPrevDepth)) >> 1);
              // - handle case diff is in base cell
              if (ddm > prevd) {
                ddm = prevd;
              }
              // - go up to common depth
              for (; ddm > 0; prevh >>= 2, ddm--, prevd--) {
                for (; (prevh & 3L) < 3L; prevh++) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
              prevh++;
              // go down!
              for (ddm = dl - prevd; ddm >= 0; prevh <<= 2, ddm--) {
                for (long targetHash = hl >>> (ddm << 1); prevh < targetHash; prevh++) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
            } else {
              // - go up to common depth
              for (; dr < prevd && (prevh & 3L) == 3L; prevh >>= 2, prevd--); // go up to first level to be completed
              for (; dr < prevd; prevh >>= 2, prevd--) { // complete going up
                for (; (++prevh & 3L) <= 3L; ) {
                  result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
                }
              }
            }
          }
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currRight = rightIt.next();  dr = currRight.getDepth();  hr = currRight.getHash();
        }
      } else { // like OR, but ignore common cells
        assert dl == dr;
        if (hl < hr) {
          result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          }
          currLeft =  leftIt.next();  dl = currLeft.getDepth();  hl = currLeft.getHash();
        } else if (hr < hl) {
          result[iResult++] = currRight.getMOCEncodedHash();
          if (!rightIt.hasNext()) {
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currRight = rightIt.next(); dr = currRight.getDepth(); hr = currRight.getHash();
        } else {
          assert hl == hr;
          // for OR, not for XOR: result[iResult++] = currLeft.getMOCEncodedHash();
          if (!leftIt.hasNext()) {
            iResult = completeWith(rightIt, result, iResult);
            break;
          } else if (!rightIt.hasNext()) { 
            iResult = completeWith(leftIt, result, iResult);
            break;
          }
          currLeft  =  leftIt.next();  dl =  currLeft.getDepth(); hl =  currLeft.getHash();
          currRight = rightIt.next();  dr = currRight.getDepth(); hr = currRight.getHash();
        }
      }
    }
    return newFrom(depthMax, result, iResult);
  }
  
  public HealpixNestedMOC not() { //  = complement
    final long[] result = new long[3 * this.size() + 12]; // Worst case: only 1 sub-cell by cell in the MOC (+11 for depth 0)
    final Iterator<CurrentCellAccessor> it = this.iterator();
    // If empty MOC 
    if (!it.hasNext()) {
      for (int i = 0; i < 12; i++) {
        result[i] = encodeHash4MOC(0, i, this.depthMax);
      }
      return new HealpixNestedMOC(this.depthMax, result);
    }
    // Else
    CurrentCellAccessor v = it.next();
    int iResult = 0, d = v.getDepth(), prevd = 0;
    long h = v.getHash(), prevh = 0;
    // Before first: go down
    for (int currd = 0, dd = d; dd >= 0; prevh <<= 2, dd--, currd++) {
      for (long targetHash = h >>> (dd << 1); prevh < targetHash; prevh++) {
        result[iResult++] = encodeHash4MOC(currd, prevh, this.depthMax);
      }
    }
    // Between first and last
    for (prevd = d, prevh = h; it.hasNext(); prevd = d, prevh = h) {
      v = it.next(); d = v.getDepth(); h = v.getHash();
      // go up (if needed)!
      // - prepare variable
      int dd = d - prevd;
      long hAtPrevDepth;
      if (dd < 0) { // case previous hash deeper that current hash
        int nBits = (-dd) << 1;
        hAtPrevDepth = h << nBits;
      } else {      // case current hash deeper that previous hash, need to go up?
        hAtPrevDepth = h >>> (dd << 1);
      }
      dd = ((63 - Long.numberOfLeadingZeros(prevh ^ hAtPrevDepth)) >> 1);
      // - handle case diff is in base cell
      if (dd > prevd) {
        dd = prevd;
      }
      // - go up to common depth
      for (; dd > 0; prevh >>= 2, dd--, prevd--) {
        for (; (prevh & 3L) < 3L; prevh++) {
          result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
        }
      }
      prevh++;
      // go down!
      for (dd = d - prevd; dd >= 0; prevh <<= 2, dd--) {
        for (long targetHash = h >>> (dd << 1); prevh < targetHash; prevh++) {
          result[iResult++] = encodeHash4MOC(prevd, prevh, this.depthMax);
        }
      }
    }
    // After last
    // - go up to depth 0
    for (; d > 0; h >>= 2, d--) {
      for (; (h & 3L) < 3L; h++) {
        result[iResult++] = encodeHash4MOC(d, h, this.depthMax);
      }
    }
    h++;
    // - complete till base cell 11
    for (; h < 12; h++) {
      result[iResult++] = encodeHash4MOC(0, h, this.depthMax);
    }
    return newFrom(this.depthMax, result, iResult);
  }
  
  // minus // = a XOR (a AND b) (like xor, but do not put cells of b)
  
  // toRanges();
  
  private static final HealpixNestedMOC newFrom(final int depthMax, final long[] a, final int i) {
    return i / (float) a.length > 0.75 ?
        new HealpixNestedMOC(depthMax, a, i)
        : new HealpixNestedMOC(depthMax, Arrays.copyOf(a, i));
  }
    
  /**
   * Create a MOC considering that the given array is already a MOC: i.e. it is sorted (ASC order)
   * and do not contains duplicate or small cells included into larger one's.
   * WARNING: the array is used internally, so it must not be modified by an external reference!
   * use {@code Arrays.copy()} is you are not sure! 
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedMOC createUnsafe(final int mocDepth, final long[] mocCells) {
    return new HealpixNestedMOC(mocDepth, mocCells);
  }

  /**
   * Same as {@link #createUnsafe(int, long[])} except that not we do not use the full array.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @param toIndex the index of the last element (exclusive) to be considered in the moc
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedMOC createUnsafe(final int mocDepth, final long[] mocCells, final int toIndex) {
    return new HealpixNestedMOC(mocDepth, mocCells);
  }

  /**
   * Same a {@link #createUnsafe(int, long[])} except that the properties (array sorted,
   * no duplicates, no cell included into an other one) is checked.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedMOC createCheck(final int mocDepth, final long[] mocCells) {
    checkIsSortedNoDuplicateNoOverlap(mocCells);
    return createUnsafe(mocDepth, mocCells);
  }

  /**
   * Same as {@link #createCheck(int, long[])} except that not we do not use the full array.
   * @param mocDepth the depth of the MOC
   * @param mocCells the array representation of the MOC
   * @param toIndex the index of the last element (exclusive) to be considered in the moc
   * @return the MOC object storing internally the array
   */
  public static HealpixNestedMOC createCheck(final int mocDepth, final long[] mocCells, final int toIndex) {
    checkIsSortedNoDuplicateNoOverlap(mocCells);
    return createUnsafe(mocDepth, mocCells, toIndex);
  }

  public static HealpixNestedBMOC create(int mocDepth, long[] mocCells) {
    // also make a version taking a sortedIterator on hash of depth < deptMax
    throw new Error("Method not yet implemented!");
    // sort 
    // remove duplicate
    // if ambiguity (flag), generate an error
  }

  @Override
  public Iterator<CurrentCellAccessor> iterator() {
    return new Iter() {
      @Override
      protected CurrentCellAccessor returnThis() {
        return this;
      }
      @Override
      protected void calledInNext(int i) { }
      @Override
      public void remove() {
       throw new UnsupportedOperationException();
      }
    };
  }

}
