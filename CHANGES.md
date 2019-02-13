
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



