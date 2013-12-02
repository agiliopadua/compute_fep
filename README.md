compute fep
===========

Free energy perturbation for LAMMPS, with soft-core pair potentials.


Installation
------------

1. Copy the `src/USER-FEP` directory into `lammps/src/USER-FEP`

2. Add `user-fep` to the `PACKUSER` variable in `lammps/src/Makefile`.

3. Compile (the `kspace` package is required):

    make yes-user-fep

    make machine


Documentation
-------------

The `compute fep`, a modified `fix adapt` and the soft-core `pair
styles` have their documentation pages in the the `doc` directory.

Some examples are provided.

