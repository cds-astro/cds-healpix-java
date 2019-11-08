
v0.25 (2019-11-08)
==================

### Added

* Add Decorator to transform a BMOC iterator into a Range iterator (to be tested with Aladin!)

### Changes

* Use a correct and more robust algo for elliptical cones

### Bug correction

* Fix test error due to elliptical cone
* Fix cone search error (not enought space allocated)


v0.24 (2019-07-15)
==================

### Bug correction

* Fix the CLI code (see closed Issue on Github)
* Fix ultra compressed BMOC

### WARNING

* Method "externalEdges" probably bugged: compare restuls with Rust!!
* BMOC logical operation are probably bugged: compare results with Rust!!


v0.23 (2019-03-21)
==================

### Changes

* Change the code of elliptical cone (still to be tested, I am wainting for an implementation in Aladin!)

### Bug correction

* Fix a trivial bug on method toRing

v0.22 (2019-03-04)
==================

### Bug correction

* Now supports longitude lower than 0 or higher than twice PI in polygone vertices

### Added

* Add support for elliptical cone


v0.21 (2019-02-28)
==================

### Bug correction

* Fix the computation of the number of rings when a cone contains a pole.
* Remove dependencies to Java 7 (nio.Files, Double.isFinite).


v0.2 (2019-02-13)
=================

### Bug correction

* Fix the exact cone solution (+ regression tests addded)
* Minimum Enclosing Cone possibly bugged for polygons larger than an hemisphere (replaced by another bounding cone).


### Added

* In polygon, we replaced the Minimum Enclosing Cone by a faster (but less accurate) bounding cone.
* Exact polygon solution (+ regression tests added)
* First version of BMOC logical operations added (but still to be debugged!!)
* Utra compact MOCs added (based on an implicit data-structure)
* First version of a CLI

### WARNINGs

* BMOC logical operations not yet debugged
* The source code is still ALL BUT CLEAN: a lot of cleaning still needed!!


v0.1 (2018-09-04)
=================

- Initial release



