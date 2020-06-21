package cds.healpix;

import java.util.Arrays;

public class MbtTest {

	/*
	public static void main( String[] args ) {
	
	Healpix.getNested( 6 )
    .newConeComputer( 0.9886484667476922 )
    .overlappingCells( 4.646003386280478, -0.41602905223485265 );

	}*/
	

  public static void main( String[] args ) {
	double[][] vertices = { 
        {4.695790732486566, 0.8889150097255261},
        {4.4682913488764395, 0.24417985540748033},
        {4.2460276551603116, -0.1307183283958091},
        {3.915276622719796, -0.42958085596528983},
        {3.159957098738856, -0.5694515052059677},
        {3.4232653287700523, -0.8411539982726408},
        {3.5624631270616547, -1.0765021986213243},
        {3.5643253154494583, -1.3247502355595913},
        {7.889934078990009, -1.4138979988201017},
        {5.093052102418385, -1.0463233540993349},
        {5.212197941830804, -0.6180885659230597},
        {5.375096075892637, -0.2088528564758103},
        {5.776366851387001, 0.33400642074068165},
        {5.4196187097295985, 0.35367158642311125},
        {5.171289457253199, 0.45115053076437905},
        {4.947552986866946, 0.6074987013677716},
     };
	
	
     Healpix.getNested( 6 ).newPolygonComputer().overlappingCells( vertices );
  }
  
  
/*public class MbtTest {
	public static void checkHash( int depth, double lon, double lat ) {
        long h0 = Healpix.getNested( depth ).hash( lon, lat );
        long h1 = Healpix.getNestedFast( depth ).hash( lon, lat );
        System.out.println( "(" + lon + "," + lat + "):\t-> " 
                          + h0 + "\t" + h1 + "\t"
                          + ( h0 == h1 ? "" : "!!!" ) );
    }
    public static void main( String[] args ) {
        int depth = 3; 
        for ( int ilon = 1; ilon >= -1; ilon-- ) {
            for ( int ilat = 1; ilat >= -1; ilat-- ) {
                checkHash( depth, 0.5 * ilon, 0.5 * ilat );
            }
        }
    }
}*/

/*

    public static void checkVertices( int depth, long hash ) {
        VerticesAndPathComputer vpc0 = Healpix.getNested( depth );
        VerticesAndPathComputer vpc1 = Healpix.getNestedFast( depth );
        CompassPoint.Cardinal east = CompassPoint.Cardinal.E;
        double[][] vs0 = vpc0.pathAlongCellEdge( hash, east, true, 1 );
        double[][] vs1 = vpc1.pathAlongCellEdge( hash, east, true, 1 );
        System.out.println( depth + "\t" + hash + ": "
                          + (Arrays.deepEquals( vs0, vs1 ) ? "ok" : "NOT OK") );
        if (!Arrays.deepEquals( vs0, vs1 )) {
        	System.out.println(depth + "\t" + hash + ": ");
        	System.out.println("Left: ");
        	print(vs0);
        	System.out.println("Right: ");
        	print(vs1);
        }
        
    }
    
    public static void print(double[][] coos) {
    	for (double[] coo : coos) {
    		System.out.println(coo[0] + ", " + coo[1]);
    	}
    }
    
    public static void main( String[] args ) {
        // for ( long h = 42236L; h < 42243L; h++ ) {
            checkVertices( 6, 42239L );
       //  }
    }
    */
}