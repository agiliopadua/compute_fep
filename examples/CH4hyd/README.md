Free Energy of Hydration of Methane
===================================

Example calculation of the free energy of hydration of methane with
LAMMPS using *compute fep* and *fix adapt*.

Methane is represented by the 5-site OPLS-AA model (1 molecule). Water
is represented by the 3-site SPC/E model (360 molecules). Interactions
of sites that are being created or deleted are treated using soft-core
verions of the Lennard-Jones and Coulomb potentials (*pair
lj/cut/coul/long/soft*) in order to avoid singularities.

The following directories contain input files and results for
calculations using free-energy perturbation (FEP) and
finite-difference thermodynamic integration (FDTI):

- `fep01` Calculation using FEP, multi-stage creation of a methane
  molecule. Results in `fep01.lmp`

- `fep10` Calculation using FEP, multi-stage deletion of a methane
  molecule. Results in `fep10.lmp`

- `fdti01` Calculation using FDTI, creation of a methane
  molecule. Results in `fdti01.lmp`

- `fdti10` Calculation using FDTI, deletion a methane
  molecule. Results in `fdti10.lmp`

The free-energy profiles can be observed by plotting the values in the
third column of the results files. The Python scripts `fep.py` and
`fdti.py` found in the `tools` directory can be used to calculate the
free-energy differences corresponding to the above transformations:

`fdti.py 300 0.002 < fdti01.lmp`

`fdti.py 300 0.002 < fdti10.lmp`

`fep.py 300 < fep01.lmp`

`fep.py 300 < fep10.lmp`

The outputs are in kcal/mol and can be compared with the experimental
value of 2.0 kcal/mol, or with a simulation value from the literature
obtained with the same force field models used here: 2.27 kcal/mol
[MR Shirts, VS Pande, J Chem Phys 122 (2005) 134508].


