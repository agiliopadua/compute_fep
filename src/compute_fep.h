/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Agilio Padua (ICCF,UBP,CNRS)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(fep,ComputeFEP)

#else

#ifndef COMPUTE_FEP_H
#define COMPUTE_FEP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeFEP : public Compute {
 public:
  int chgflag; 

  ComputeFEP(class LAMMPS *, int, char **);
  ~ComputeFEP();
  void init();
  double compute_scalar();
  void compute_vector();

 private:
  int npert;
  int anypair;
  int tailflag;
  double temp_fep;

  double **f_orig;
  double eng_vdwl_orig,eng_coul_orig;
  double pvirial_orig[6];
  double *peatom_orig,**pvatom_orig;
  double energy_orig;
  double kvirial_orig[6];
  double *keatom_orig,**kvatom_orig;

  struct Perturb {
    int which;
    char *pstyle,*pparam;
    int ilo,ihi,jlo,jhi;
    int pdim;
    double delta;
    double *q_orig;
    double **array,**array_orig;
    int aparam;
  };

  Perturb *perturb;

  class Pair *pair;

  void change_params();
  double compute_epair();
  void restore_params();
  void backup_accumulators();
  void restore_accumulators();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for compute fep does not exist

Self-explanatory.

E: Variable for compute fep is invalid style

Self-explanatory.

E: Compute fep pair style does not exist

Self-explanatory.

E: Energy was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
