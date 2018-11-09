
About
-----

CDS implementation in Java of the HEALPix tesselation.

License
-------

This software is released under the 3-Clause BSD license.
See the license [file](LICENSE.md).  


Install
-------

You need [ant](https://ant.apache.org/) to be installed on your machine to test and build the project.

Then, in the root of the project (the directory containing the build.xml file, type
```bash
> ant -p
```
to list possible actions.  
You should get:
```
CDS Healpix library build.xml
Main targets:

 make.all    Compile all the software, run tests, and create the .jar
 make.build  Compile main classes, to test if the soft compiles
 make.doc    Compile all the software and create the .jar
 make.jar    Compile main classes and create both the jar and the fatjar
 make.test   Compile all classes (main an test), and run tests
Default target: print_properties
```

To build the jar, simply type
```bash
>  ant make.jar
```

To build the javadoc:
```bash
>  ant make.doc
```
Then open it in you favourite browser, e.g.:
```bash
>  firefox target/docs/apidocs/index.html
```

The test contains performances tests, run all of them with:
```bash
>  ant make.test
```

Usage examples
--------------

#### Cell number from equatorial coordinates at a given depth (order) 

```java
// Inputs
int depth = 8;
double ra  = Math.toRadians(10.509734);
double dec = Math.toRadians(21.657381);

// Result
HealpixNested hn = Healpix.getNested(depth);
long cellNumber = hn.hash(ra, dec);

// Or
HealpixNestedFast hnf = Healpix.getNestedFast(depth);
cellNumber = hnf.hash(ra, dec);

// In a multi-threaded environment, use one HashComputer per thread
HashComputer hc = hn.newHashComputer();
cellNumber = hc.hash(ra, dec);
```

#### Coordinates on the sky of the center (in the projection plane) of a cell

```java
// Static imports
import static cds.healpix.VerticesAndPathComputer.LON_INDEX;
import static cds.healpix.VerticesAndPathComputer.LAT_INDEX;

// Inputs
int depth = 8;
long cellNumber = 12394L;

// Result
HealpixNested hn = Healpix.getNested(depth);
double[] centerCoos = hn.center(cellNumber);
double raDeg  = Math.toDegrees(centerCoos[LON_INDEX]);
double decDeg = Math.toDegrees(centerCoos[LAT_INDEX]);

// Again, in a multi-threaded environment, use one VerticesAndPathComputer per thread
VerticesAndPathComputer vpc = hn.newVerticesAndPathComputer();
centerCoos = vpc.center(cellNumber);
raDeg  = Math.toDegrees(centerCoos[LON_INDEX]);
decDeg = Math.toDegrees(centerCoos[LAT_INDEX]);
```

#### Coordinates on the sky of the 4 vertices of a cell

```java
// Static import
import static cds.healpix.VerticesAndPathComputer.ALL_CARDINAL_POINTS;

// Inputs
int depth = 8;
long cellNumber = 12394L;

// Result
HealpixNested hn = Healpix.getNested(depth);
EnumMap<Cardinal, double[]> vertices = vertices(cellNumber, ALL_CARDINAL_POINTS);

// and then
for (double[] vertexCoos : vertices) {
    ...
}

// or
double[] northVertexCoos = vertices.get(Cardinal.N);
double[] eastVertexCoos  = vertices.get(Cardinal.E);
double[] southVertexCoos = vertices.get(Cardinal.S);
double[] westVertexCoos  = vertices.get(Cardinal.W);


// If you are only interested in a single vertex
northVertexCoos = hn.vertex(Cardinal.N);

// If you are only interested in East and West vertices
vertices = vertices(cellNumber, EnumSet.of(Cardinal.E, Cardinal.W));
```

#### Cell neighbours

```java
// Inputs
int depth = 8;
long cellNumber = 12394L;

// Results
HealpixNested hn = Healpix.getNested(depth);
NeighbourList neigList = hn.neighbours(cellNumber);

// Iterates
FlatHashIterator hIt = neigList.iterator();
while (hIt.hasNext()) {
    System.out.println("Neighbour cell number: " + hIt.next());
}

// Or access by MainWinds
long northeastCellNumber = neigList.get(MainWind.NE);
```

#### Cells in a cone

```java
// Inputs
int depth = 8;
double coneCenterRa  = Math.toRadians(10.509734);
double coneCenterDec = Math.toRadians(21.657381);
double coneRadiusRad = Math.toRadians(1.0);

// Resulting BMOC (MOC + flag telling is a cell is fully or partially overlapped by the cone)
HealpixNested hn = Healpix.getNested(depth);
HealpixNestedFixedRadiusConeComputer cc = hn.newConeComputer(coneRadiusRad);
HealpixNestedBMOC bmoc = cc.overlappingCells(coneCenterLonRad, coneCenterLatRad);

// BMoc view
for (HealpixNestedBMOC.CurrentValueAccessor cell : bmoc) {
    System.out.println("Cell number: " + cell.getHash() 
        + "; cell depth:" + cell.getDepth()
        + "; fully covered: " + cell.isFull());
}

// Flat view at depth 8:
FlatHashIterator hIt = bmoc.flatHashIterator();
while (hIt.hasNext()) {
    System.out.println("Cell number at depth 8: " + hIt.next());
}
```

#### Cone search for cross-matches

```java
// Inputs
int depth = 8;
double xmatchRadius = Math.toRadians(1.0 / 3600.0); // 1 arcsec

// Results
HealpixNested hn = Healpix.getNested(depth);
HealpixNestedFixedRadiusCone4XMatch xm = hn.newConeComputer4Xmatch(xmatchRadius);

long[] cellsNum = new long[9];
int indexMax = 0;
for (double[] coosInRadians : listOfPosToXmatch) {
  indexMax = xm.overlappingCells(coosInRadians[0], coosInRadians[1], cellsNum);
  for (int ic = 0; i < indexMax; i++) {
    // look at cellsNum[ic], ...
  }
}
```




