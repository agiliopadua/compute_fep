compute fep
===========

_[Agilio Padua](http://tim.univ-bpclermont.fr/apadua)_

Free energy perturbation for LAMMPS, with soft-core pair potentials.


Installation
------------

The `compute fep` is installed as a user-package.

1. Copy the `src/USER-FEP` directory to `lammps/src/USER-FEP`.

2. You may want to keep a backup of the original `fix_adapt.h` and
   `fix_adapt.cpp`, although the version included here is more flexible.

2. Add `user-fep` to the `PACKUSER` variable in `lammps/src/Makefile`

3. Compile (the `kspace` package is required)

        make yes-user-fep

        make 'machine'


Documentation
-------------

The `compute fep`, a modified `fix adapt` and the soft-core `pair
styles` have their documentation pages in the the `doc` directory.

Some examples are provided.
