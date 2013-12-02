compute_fep
===========

Free energy perturbation for LAMMPS, with soft-core pair potentials.


Installation
------------

1. Copy the `src/USER-FEP` directory into `lammps/src/USER-FEP`

2. Add `user-fep` to the `PACKUSER` variable in `lammps/src/Makefile`.

3. Compile (the `molecule` and `kspace` packages are required):

    make yes-user-fep

    make machine


Documentation
-------------

The `compute fep`, the modified `fix adapt` and the soft-core `pair
styles` have documentation pages in the the `doc` directory.

Some examples of usage are provided.

