package cds.healpix;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.EnumMap;
import java.util.Map;
import java.util.Scanner;

import cds.healpix.CompassPoint.Cardinal;

import static cds.healpix.VerticesAndPathComputer.LON_INDEX;
import static cds.healpix.VerticesAndPathComputer.LAT_INDEX;
import static cds.healpix.VerticesAndPathComputer.ALL_CARDINAL_POINTS;

/**
 * Simple and Basic CLI (Command Line Interface) for quick testing.
 * 
 * @author F.-X. Pineau
 *
 */
public class HealpixCLI {


  private static final String usage() {
    final StringBuilder s = new StringBuilder();
    s.append("HealpixCLI action params").append('\n');
    s.append("ACTION").append('\n');
    s.append("  help, -h, -help, --help      print this message").append('\n');
    s.append("  quit, exit                   exits").append('\n');
    s.append("  hash     DEPTH RA DEC        compute a cell number at the given depth from the given coords").append('\n');
    s.append("  center   DEPTH HASH          compute the coordinates (in degrees) of a cell center").append('\n');
    s.append("  vertices DEPTH HASH          compute the coordinates (in degrees) of the 4 cell vertices").append('\n');
    s.append("  neigh    DEPTH HASH          provides the cell number of each neighbour of the given cell").append('\n');
    s.append("  cone     DEPTH RA DEC RADIUS provides the list of cells overlaped by a cone").append('\n');
    s.append("PARAMS").append('\n');
    s.append("  DEPTH    healpix layer depth (in [0, 29]").append('\n');
    s.append("  RA       right ascension, in [-180, 360] degrees").append('\n');
    s.append("  DEC      declination, in [-90, 90] degrees").append('\n');
    s.append("  HASH     an HEALPix hash value, i.e. a cell number").append('\n');
    s.append("  RADIUS   a radius, in degrees").append('\n');
    return s.toString();
  }

  private static final void printUsage() {
    System.out.println(usage());
  }

  private static final Deque<String> args2stack(final String... args) {
    final Deque<String> stack = new ArrayDeque<String>(args.length);
    for (final String arg : args) {
      stack.addLast(arg);
    }
    return stack;
  }

  private static int popDepth(final Deque<String> stack) {
    if (stack.isEmpty()) {
      throw new IllegalArgumentException("DEPTH parameter is missing");
    }
    final String depthStr = stack.pop();
    int depth = 0;
    try {
      depth = Integer.parseInt(depthStr);
    } catch(NumberFormatException e) {
      throw new IllegalArgumentException("DEPTH parameter. Expected: integer; Actual: " + depthStr);
    }
    return depth;
  }

  private static long popHash(final Deque<String> stack) {
    if (stack.isEmpty()) {
      throw new IllegalArgumentException("HASH parameter is missing");
    }
    final String hStr = stack.pop();
    long h = 0;
    try {
      h = Long.parseLong(hStr);
    } catch(NumberFormatException e) {
      throw new IllegalArgumentException("DEPTH parameter. Expected: integer; Actual: " + hStr);
    }
    return h;
  }

  private static double popDouble(final Deque<String> stack, String paramName) {
    if (stack.isEmpty()) {
      throw new IllegalArgumentException(paramName + " parameter is missing");
    }
    final String dStr = stack.pop();
    double d = 0;
    try {
      d = Double.parseDouble(dStr);
    } catch(NumberFormatException e) {
      throw new IllegalArgumentException(paramName + " parameter. Expected: double; Actual: " + dStr);
    }
    return d;
  }

  private static double popRA(final Deque<String> stack) {
    return Math.toRadians(popDouble(stack, "RA"));
  }

  private static double popDec(final Deque<String> stack) {
    return Math.toRadians(popDouble(stack, "DEC"));
  }

  private static double popRadius(final Deque<String> stack) {
    return Math.toRadians(popDouble(stack, "RADIUS"));
  }

  private static final void execCatchErr(final Deque<String> stack) {
    try {
      exec(stack);
    } catch(Exception e) {
      System.out.println("Error: " + e.getMessage());
    } catch(Error e) {
      System.out.println("Error: " + e.getMessage());
    }
  }

  /* REMOVED TO STAY COMPATIBLE WITH Java6
   private static final void exec(final Deque<String> stack) {
    final String action = stack.pop();
    switch (action) {
    case "help": case "-h": case "-help": case "--help":
      printUsage();
      break;
    case "quit": case "exit":
      System.exit(0);
    case "hash": {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      long cellNumber = hn.hash(popRA(stack), popDec(stack));
      System.out.println(cellNumber);
      break;
    }
    case "center": {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      double[] centerCoos = hn.center(popHash(stack));
      double raDeg  = Math.toDegrees(centerCoos[LON_INDEX]);
      double decDeg = Math.toDegrees(centerCoos[LAT_INDEX]);
      System.out.println(raDeg + " " + decDeg);
      break;
    }
    case "vertices": {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      EnumMap<Cardinal, double[]> vertices = hn.vertices(popHash(stack), ALL_CARDINAL_POINTS);
      for (final Map.Entry<Cardinal, double[]> e : vertices.entrySet()) {
        final Cardinal c = e.getKey();
        final double[] coos = e.getValue();
        double raDeg  = Math.toDegrees(coos[LON_INDEX]);
        double decDeg = Math.toDegrees(coos[LAT_INDEX]);
        System.out.println(c + ": " + raDeg + " " + decDeg);
      }
      break;
    }
    case "neigh": {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      NeighbourList neigList = hn.neighbours(popHash(stack));
      for (int i = 0; i < neigList.size(); i++) {
        System.out.println(neigList.getDirection(i) + ": " + neigList.get(i));
      }
      break;
    }
    case "cone": {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      double coneCenterRa  = popRA(stack);
      double coneCenterDec = popDec(stack);
      HealpixNestedFixedRadiusConeComputer cc = hn.newConeComputer(popRadius(stack));
      HealpixNestedBMOC bmoc = cc.overlappingCells(coneCenterRa, coneCenterDec);
      for (HealpixNestedBMOC.CurrentValueAccessor cell : bmoc) {
        String isFull = cell.isFull() ? "(f)" : "(p)";
        System.out.println(cell.getDepth()+ "/" + cell.getHash() + isFull);
      }
      break;
    }
    case "":
      break;
    default:
      throw new IllegalArgumentException("Unknown action \"" + action + "\"");
    }
  }*/
  
  private static final void exec(final Deque<String> stack) {
    final String action = stack.pop();
    System.out.println("action: " + action);
    if (   action.equals("help")  || action.equals("-h") 
        || action.equals("-help") || action.equals("--help")) {
      printUsage();
    } else if (action.equals("quit") || action.equals("exit")) {
      System.exit(0);
    } else if (action.equals("hash")) {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      long cellNumber = hn.hash(popRA(stack), popDec(stack));
      System.out.println(cellNumber);
    } else if (action.equals("center")) {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      double[] centerCoos = hn.center(popHash(stack));
      double raDeg  = Math.toDegrees(centerCoos[LON_INDEX]);
      double decDeg = Math.toDegrees(centerCoos[LAT_INDEX]);
      System.out.println(raDeg + " " + decDeg);
    } else if (action.equals("vertices")) {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      EnumMap<Cardinal, double[]> vertices = hn.vertices(popHash(stack), ALL_CARDINAL_POINTS);
      for (final Map.Entry<Cardinal, double[]> e : vertices.entrySet()) {
        final Cardinal c = e.getKey();
        final double[] coos = e.getValue();
        double raDeg  = Math.toDegrees(coos[LON_INDEX]);
        double decDeg = Math.toDegrees(coos[LAT_INDEX]);
        System.out.println(c + ": " + raDeg + " " + decDeg);
      }
    } else if (action.equals("neigh")) {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      NeighbourList neigList = hn.neighbours(popHash(stack));
      for (int i = 0; i < neigList.size(); i++) {
        System.out.println(neigList.getDirection(i) + ": " + neigList.get(i));
      }
    } else if (action.equals("cone")) {
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      double coneCenterRa  = popRA(stack);
      double coneCenterDec = popDec(stack);
      HealpixNestedFixedRadiusConeComputer cc = hn.newConeComputer(popRadius(stack));
      HealpixNestedBMOC bmoc = cc.overlappingCells(coneCenterRa, coneCenterDec);
      for (HealpixNestedBMOC.CurrentValueAccessor cell : bmoc) {
        String isFull = cell.isFull() ? "(f)" : "(p)";
        System.out.println(cell.getDepth()+ "/" + cell.getHash() + isFull);
      }
    } else if (action.equals("")) {
    } else {
      throw new IllegalArgumentException("Unknown action \"" + action + "\"");
    }
  }

  public static final void main(String[] args) {
    if (args == null || args.length == 0) { // Interactive mode
      printUsage();
      Scanner in = new Scanner(System.in);
      do {
        args = in.nextLine().split("\\s+");
        execCatchErr(args2stack(args));
      } while(true);
      // in.close(); // unreachable, we don't care since brutal exit
    } else { // One shot mode
      exec(args2stack(args));
    }
  }

}
