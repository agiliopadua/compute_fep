compute fep
===========

_[Agilio Padua](http://tim.univ-bpclermont.fr/apadua)_

Free energy perturbation for LAMMPS, with soft-core pair potentials.


Installation
------------

The `compute fep`, the `fix adapt/fep` and several soft-core `pair styles`
are included in the LAMMPS distribution (main code and also the
LAMMPS-ICMS version) as the USER-FEP package. The `kspace` package is
required.

To compile the USER-FEP package, cd to `lammps/src`:

        make yes-user-fep

        make 'machine'


Documentation
-------------

The `compute fep`, the `fix adapt/fep` and the soft-core `pair
styles` have their documentation pages in the the `doc` directory.

Some examples are provided.
