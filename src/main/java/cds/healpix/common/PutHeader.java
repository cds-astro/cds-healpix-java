// Copyright 2017-2018 - Université de Strasbourg/CNRS
// The CDS HEALPix library is developped by the Centre de Données
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

package cds.healpix.common;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.StringTokenizer;

/**
 * 
 * @author P. Fernique, modified by F.-X. Pineau
 *
 */
class PutHeader {
  String PATH = "/home/pineau/Eclipse/Healpix/src";
  String [] EXCEPTION = { "savot","astrores","table" };
  String HEAD =
          "// Copyright 2017-2018 - Université de Strasbourg/CNRS\n" +
          "// The CDS HEALPix library is developped by the Centre de Données\n" +
          "// astronomiques de Strasbourgs (CDS) from the following external papers:\n" +
          "//  - [Gorsky2005]     - \"HEALPix: A Framework for High-Resolution Discretization and\n" +
          "//                       Fast Analysis of Data Distributed on the Sphere\"\n" +
          "//                       http://adsabs.harvard.edu/abs/2005ApJ...622..759G\n" +
          "//  - [Calabretta2004] - \"Mapping on the HEALPix grid\"\n" +
          "//                       http://adsabs.harvard.edu/abs/2004astro.ph.12607C\n" +
          "//  - [Calabretta2007] - \"Mapping on the HEALPix grid\"\n" +
          "//                       http://adsabs.harvard.edu/abs/2007MNRAS.381..865C\n" +
          "//  - [Reinecke2015]   - \"Efficient data structures for masks on 2D grids\"\n" +
          "//                       http://adsabs.harvard.edu/abs/2015A&A...580A.132R\n" +
          "// It is distributed under the terms of the BSD License 2.0\n" +
          "//\n" +
          "// This file is part of the CDS HEALPix library.\n" +
          "//\n";

  PutHeader() {
    int n = putHead(new File(PATH), 0);
    System.out.println("Entête maj dans "+n+" fichier"+(n>1?"s":""));
  }

  int putHead(File dir, int depth) {
    int n=0;
    File [] files = dir.listFiles();
    for( File f: files ) {
      try {
        String name = f.getCanonicalPath();
        if( f.isDirectory() ) {
          // Peut être une branche à ne pas prendre en compte ?
          if( depth==0 ) {
            int off = name.lastIndexOf('/');
            if( off==-1 ) off = name.lastIndexOf('\\');
            String d = name.substring(off+1);
            boolean trouve=false;
            for( String ex : EXCEPTION ) {
              if( d.equals(ex) ) { trouve=true; break; }
            }
            if( trouve ) continue;
          }

          n += putHead(f,depth+1);
        }

        if( !f.getCanonicalPath().endsWith(".java") ) continue;
        insertHead(f);
        n++;
        //            if( true ) {
        //               System.out.println("C'est fait pour "+f.getCanonicalPath());
        //               if( n==2  ) System.exit(1);
        //            }
      } catch( Exception e ) { e.printStackTrace(); }
    }
    return n;
  }

  void insertHead(File f) throws Exception {
    //      long time = f.lastModified();
    String tmp = f.getCanonicalPath()+".tmp";
    File g = new File(tmp);
    BufferedWriter t = null;
    BufferedReader s = null;
    try {
      s = new BufferedReader( new InputStreamReader( new FileInputStream(f)));
      if( g.exists() ) g.delete();
      t = new BufferedWriter( new OutputStreamWriter( new FileOutputStream(g)));
      writeHead(t);
      skipPreviousHead(s);
      writeContent(s,t);
      t.close();
      s.close();

      f.delete();
      //         g.setLastModified(time);
      g.renameTo(f);
    }
    catch( Exception e) {
      System.out.println("Problème sur "+f.getCanonicalPath()+" ["+e.getMessage()+"] -> ignoré");
    }
    finally {
      if( t!=null ) t.close();
      if( s!=null ) s.close();
    }
  }

  void writeContent(BufferedReader s,BufferedWriter t) throws Exception {
    String s1;
    while( (s1=s.readLine())!=null ) {
      t.write(s1);
      t.newLine();
    }
  }

  void skipPreviousHead(BufferedReader s) throws Exception {
    while( true ) {
      s.mark(2048);
      String s1 = s.readLine();
      if( s1==null ) break;
      if( s1.trim().length()>0 && !s1.startsWith("//") ) {
        s.reset();
        return;
      }
      s1 = s1.trim();
    }
  }

  void writeHead(BufferedWriter t) throws Exception {
    StringTokenizer st = new StringTokenizer(HEAD, "\n");
    while( st.hasMoreTokens() ) {
      t.write(st.nextToken());
      t.newLine();
    }
    t.newLine();
  }

  public static void main(String[] args) {
    try {
      new PutHeader();
    } catch( Exception e) { e.printStackTrace(); }

  }

} 
