
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

