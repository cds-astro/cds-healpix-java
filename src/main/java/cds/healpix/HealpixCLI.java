package cds.healpix;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Scanner;

public class HealpixCLI {

  
  private static final String usage() {
    final StringBuilder s = new StringBuilder();
    s.append("HealpixCLI action params").append('\n');
    s.append("ACTIONS").append('\n');
    s.append("  help, -h, -help, --help    print this message").append('\n');
    s.append("  quit, exit                 exits").append('\n');
    s.append("  hash DEPTH RA_DEG DEC_DEG  compute a cell number at the given depth from the given coords").append('\n');
    return s.toString();
  }
  
  private static final void printUsage() {
    System.out.println(usage());
  }
  
  /*private static final String[] getUserInput() {
    Scanner in = new Scanner(System.in);
    final String[] args = in.nextLine().split("\\s+");
    in.close();
    return args;
  }*/
  
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
  
  private static final void execCatchErr(final Deque<String> stack) {
    try {
      exec(stack);
    } catch(Exception e) {
      System.out.println("Error: " + e.getMessage());
    } catch(Error e) {
      System.out.println("Error: " + e.getMessage());
    }
  }
  
  private static final void exec(final Deque<String> stack) {
    final String action = stack.pop();
    switch (action) {
    case "help": case "-h": case "-help": case "--help":
      printUsage();
      break;
    case "quit": case "exit":
      System.exit(0);
    case "hash":
      HealpixNested hn = Healpix.getNested(popDepth(stack));
      long cellNumber = hn.hash(popRA(stack), popDec(stack));
      System.out.println(cellNumber);
      break;
    default:
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
    } else { // One shot mode
      exec(args2stack(args));
    }
  }
  
}
