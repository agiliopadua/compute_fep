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
   Contributing author: Agilio Padua (ICCF UBP CNRS)
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "mpi.h"
#include "comm.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "compute_fep.h"

using namespace LAMMPS_NS;

enum{PAIR,ATOM};
enum{CHARGE};

#define FEP_DEBUG

/* ---------------------------------------------------------------------- */

ComputeFEP::ComputeFEP(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal number of arguments in compute fep");

  scalar_flag = 1;

  temp_fep = force->numeric(FLERR,arg[3]);

  // count # of perturbations

  npert = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal pair attribute in compute fep");
      npert++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal atom attribute in compute fep");
      npert++;
      iarg += 4;
    } else break;
  }

  if (npert == 0) error->all(FLERR,"Illegal syntax in compute fep");
  perturb = new Perturb[npert];

  // parse keywords

  npert = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      perturb[npert].which = PAIR;
      int n = strlen(arg[iarg+1]) + 1;
      perturb[npert].pstyle = new char[n];
      strcpy(perturb[npert].pstyle,arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      perturb[npert].pparam = new char[n];
      strcpy(perturb[npert].pparam,arg[iarg+2]);
      force->bounds(arg[iarg+3],atom->ntypes,
                    perturb[npert].ilo,perturb[npert].ihi);
      force->bounds(arg[iarg+4],atom->ntypes,
                    perturb[npert].jlo,perturb[npert].jhi);
      perturb[npert].delta = force->numeric(FLERR,arg[iarg+5]);
      npert++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      perturb[npert].which = ATOM;
      if (strcmp(arg[iarg+1],"charge") == 0) {
        perturb[npert].aparam = CHARGE; 
        force->bounds(arg[iarg+2],atom->ntypes,
                      perturb[npert].ilo,perturb[npert].ihi);
      } else error->all(FLERR,"Illegal atom argument in compute fep");
      perturb[npert].delta = force->numeric(FLERR,arg[iarg+3]);
      npert++;
      iarg += 4;
    } else break;
  }

  // optional keywords

  tailflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"tail") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal optional keyword in compute fep");
      if (strcmp(arg[iarg+1],"no") == 0) tailflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) tailflag = 1;
      else error->all(FLERR,"Illegal optional keyword in compute fep");
      iarg += 2;
    } else error->all(FLERR,"Illegal optional keyword in compute fep");
  }

#ifdef FEP_DEBUG
  if (comm->me == 0 && screen) {
    fprintf(screen, "###FEP temperature = %7.2f K\n", temp_fep);
    fprintf(screen, "###FEP tail %s\n", (tailflag ? "yes":"no"));
    for (int m = 0; m < npert; m++) {
      Perturb *pert = &perturb[m];
      if (pert->which == PAIR)
        fprintf(screen, "###FEP %s %s %d-%d %d-%d %f\n", pert->pstyle, pert->pparam,
                pert->ilo, pert->ihi, pert->jlo, pert->jhi, pert->delta);
      else if (pert->which == ATOM)
        fprintf(screen, "###FEP %d-%d charge %f\n", pert->ilo, pert->ihi, pert->delta);
    }
  }
#endif

  // allocate pair style arrays

  int ntype = atom->ntypes;
  int natom = atom->nlocal + atom->nghost;
  for (int m = 0; m < npert; m++) {
    if (perturb[m].which == PAIR)
      memory->create(perturb[m].array_orig,ntype+1,ntype+1,"fep:array_orig");
    else if (perturb[m].which == ATOM) {
      memory->create(perturb[m].q_orig,natom+1,"fep:q_orig");
    }
  }

  // allocate arrays for force, energy, virial backups

  natom = atom->natoms;
  memory->create(f_orig,natom+1,3,"fep:f_orig");
  memory->create(peatom_orig,natom+1,"fep:peatom_orig");
  memory->create(pvatom_orig,natom+1,6,"fep:pvatom_orig");
  if (force->kspace) {
    memory->create(keatom_orig,natom+1,"fep:keatom_orig");
    memory->create(kvatom_orig,natom+1,6,"fep:kvatom_orig");
  }

}

/* ---------------------------------------------------------------------- */

ComputeFEP::~ComputeFEP()
{
  for (int m = 0; m < npert; m++) {
    if (perturb[m].which == PAIR) {
      delete [] perturb[m].pstyle;
      delete [] perturb[m].pparam;
      memory->destroy(perturb[m].array_orig);
    } else if (perturb[m].which == ATOM)
       memory->destroy(perturb[m].q_orig);
  }
  delete [] perturb;

  memory->destroy(f_orig);
  memory->destroy(peatom_orig);
  memory->destroy(pvatom_orig);
  if (force->kspace) {
    memory->destroy(keatom_orig);
    memory->destroy(kvatom_orig);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFEP::init()
{
  int i,j;

  // setup and error checks

  anypair = 0;

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];

    if (force->pair == NULL)
      error->all(FLERR,"compute fep pair requires pair interactions");

    if (pert->which == PAIR) {

      anypair = 1;

      Pair *pair = force->pair_match(pert->pstyle,1);
      if (pair == NULL) error->all(FLERR,"compute fep pair style does not exist");
      void *ptr = pair->extract(pert->pparam,pert->pdim);
      if (ptr == NULL) 
        error->all(FLERR,"compute fep pair style param not supported");

      pert->array = (double **) ptr;

      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

      if ((strcmp(force->pair_style,"hybrid") == 0 ||
           strcmp(force->pair_style,"hybrid/overlay") == 0)) {
        PairHybrid *pair = (PairHybrid *) force->pair;
        for (i = pert->ilo; i <= pert->ihi; i++)
          for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
            if (!pair->check_ijtype(i,j,pert->pstyle))
              error->all(FLERR,"compute fep type pair range is not valid for "
                         "pair hybrid sub-style");
      }

    } else if (pert->which == ATOM) {
      if (pert->aparam == CHARGE) {
        if (!atom->q_flag)
          error->all(FLERR,"compute fep requires atom attribute charge");
      }
    }
  }

  if (tailflag) {
    if (force->pair->tail_flag == 0)
      error->all(FLERR,"Compute fep tail when pair style does not "
                 "compute tail corrections");
  }

  // make copy of original pair,atom array values

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];
    if (pert->which == PAIR) {
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
          pert->array_orig[i][j] = pert->array[i][j];
    } else if (pert->which == ATOM) {
        int *atype = atom->type;
        double *q = atom->q; 
        int *mask = atom->mask;
        int natom = atom->nlocal + atom->nghost;
        for (i = 0; i < natom; i++)
          if (atype[i] >= pert->ilo && atype[i] <= pert->ihi)
            if (mask[i] & groupbit)
              pert->q_orig[i] = q[i];         
    }
  }

#ifdef FEP_DEBUG
  // need an occasional half neighbor list
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
#endif

}

/* ---------------------------------------------------------------------- */

void ComputeFEP::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeFEP::compute_scalar()
{
  double pe0,pe1;

#ifdef FEP_DEBUG
  neighbor->build_one(list->index);
  double pe0_test = energy_pair();
#endif

  pe0 = compute_epair();

  change_params();

#ifdef FEP_DEBUG
  double pe1_test = energy_pair();
#endif

  timer->stamp();
  if (force->pair && force->pair->compute_flag) {
    force->pair->compute(1,0);
    timer->stamp(TIME_PAIR);
  }

  if (force->kspace && force->kspace->compute_flag) {
    force->kspace->compute(1,0);
    timer->stamp(TIME_KSPACE);
  }
  pe1 = compute_epair();

  restore_params();

#ifdef FEP_DEBUG
  if (comm->me == 0 && screen)
    fprintf(screen, "###FEP pe0 = %e(%e)  pe1 = %e(%e)\n",
            pe0,pe0_test,pe1,pe1_test);
#endif

  scalar = exp(-(pe1-pe0)/(force->boltz*temp_fep));
  return scalar;
}


/* ----------------------------------------------------------------------
   obtain pair energy from lammps accumulators
------------------------------------------------------------------------- */

double ComputeFEP::compute_epair()
{
  double eng, eng_pair;

  eng = 0.0;
  if (force->pair)
    eng = force->pair->eng_vdwl + force->pair->eng_coul;
  MPI_Allreduce(&eng,&eng_pair,1,MPI_DOUBLE,MPI_SUM,world);

  if (tailflag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    eng_pair += force->pair->etail / volume;
  }

  if (force->kspace) eng_pair += force->kspace->energy;

  return eng_pair;
}

/* ----------------------------------------------------------------------
   change pair,kspace,atom parameters based on variable evaluation
------------------------------------------------------------------------- */

void ComputeFEP::change_params()
{
  int i,j;

  // make copy of original force, energy, virial array values

  backup_accumulators();

  // apply perturbation to interaction parameters

  for (int m = 0; m < npert; m++) {
    Perturb *pert = &perturb[m];

    if (pert->which == PAIR) {      // modify pair parameters
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
          pert->array[i][j] = pert->array_orig[i][j] + pert->delta;
      
#ifdef FEP_DEBUG
      if (comm->me == 0 && screen) {
        fprintf(screen, "###FEP change %s %s\n", pert->pstyle, pert->pparam);
        fprintf(screen, "###FEP  I  J   old_param new_param\n");
        for (i = pert->ilo; i <= pert->ihi; i++)
          for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
            fprintf(screen, "###FEP %2d %2d %9.5f %9.5f\n", i, j, 
                    pert->array_orig[i][j], pert->array[i][j]);
      }
#endif

    } else if (pert->which == ATOM) {

      if (pert->aparam == CHARGE) {      // modify charges
        int *atype = atom->type;
        double *q = atom->q; 
        int *mask = atom->mask;
        int natom = atom->nlocal + atom->nghost;
        for (i = 0; i < natom; i++)
          if (atype[i] >= pert->ilo && atype[i] <= pert->ihi)
            if (mask[i] & groupbit)
              q[i] += pert->delta; 

#ifdef FEP_DEBUG
        if (comm->me == 0 && screen) {
          fprintf(screen, "###FEP change charge\n");
          fprintf(screen, "###FEP  atom  I   old_q     new_q\n");
          for (i = 0; i < atom->nlocal; i++)
            if (atype[i] >= pert->ilo && atype[i] <= pert->ihi)
              if (mask[i] & groupbit)
                fprintf(screen, "###FEP %5d %2d %9.5f %9.5f\n", i, atype[i],
                        pert->q_orig[i], q[i]);
        }
#endif
      }

    }
  }

  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections

  if (anypair) force->pair->reinit();
}

/* ----------------------------------------------------------------------
   restore pair,atom parameters to original values
------------------------------------------------------------------------- */

void ComputeFEP::restore_params()
{
  int i,j;

  for (int m = 0; m < npert; m++) {

    // restore pair parameters

    Perturb *pert = &perturb[m];
    if (pert->which == PAIR) {
      for (i = pert->ilo; i <= pert->ihi; i++)
        for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
          pert->array[i][j] = pert->array_orig[i][j];
#ifdef FEP_DEBUG
      if (comm->me == 0 && screen) {
        fprintf(screen, "###FEP restore %s %s\n", pert->pstyle, pert->pparam);
        fprintf(screen, "###FEP  I  J   param\n");
        for (i = pert->ilo; i <= pert->ihi; i++)
          for (j = MAX(pert->jlo,i); j <= pert->jhi; j++)
            fprintf(screen, "###FEP %2d %2d %9.5f\n", i, j, pert->array[i][j]);
      }
#endif

    } else if (pert->which == ATOM) {

      // restore charges

      int *atype = atom->type;
      double *q = atom->q; 
      int *mask = atom->mask;
      int natom = atom->nlocal + atom->nghost;
      for (i = 0; i < natom; i++)
        if (atype[i] >= pert->ilo && atype[i] <= pert->ihi)
          if (mask[i] & groupbit)
            q[i] = pert->q_orig[i]; 
      
#ifdef FEP_DEBUG
      if (comm->me == 0 && screen) {
        fprintf(screen, "###FEP restore charge\n");
        fprintf(screen, "###FEP  atom  I   q\n");
        for (i = 0; i < atom->nlocal; i++)
          if (atype[i] >= pert->ilo && atype[i] <= pert->ihi)
            if (mask[i] & groupbit)
              fprintf(screen, "###FEP %5d %2d %9.5f\n", i, atype[i], q[i]);
      }
#endif
    }
  }

  // restore force, energy, virial array values

  restore_accumulators();

  if (anypair) force->pair->reinit();
}


/* ----------------------------------------------------------------------
   backup arrays with force, energy, virial accumulators
------------------------------------------------------------------------- */

void ComputeFEP::backup_accumulators()
{
  int i;

  int natom = atom->nlocal;
  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f_orig[i][0] = f[i][0]; 
    f_orig[i][1] = f[i][1]; 
    f_orig[i][2] = f[i][2]; 
  }

  eng_vdwl_orig = force->pair->eng_vdwl;
  eng_coul_orig = force->pair->eng_coul;

  pvirial_orig[0] = force->pair->virial[0];
  pvirial_orig[1] = force->pair->virial[1];
  pvirial_orig[2] = force->pair->virial[2];
  pvirial_orig[3] = force->pair->virial[3];
  pvirial_orig[4] = force->pair->virial[4];
  pvirial_orig[5] = force->pair->virial[5];

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++)
      peatom_orig[i] = peatom[i];
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++) {
      pvatom_orig[i][0] = pvatom[i][0];
      pvatom_orig[i][1] = pvatom[i][1];
      pvatom_orig[i][2] = pvatom[i][2];
      pvatom_orig[i][3] = pvatom[i][3];
      pvatom_orig[i][4] = pvatom[i][4];
      pvatom_orig[i][5] = pvatom[i][5];
    }
  }

  if (force->kspace) {
    energy_orig = force->kspace->energy;
    kvirial_orig[0] = force->kspace->virial[0];
    kvirial_orig[1] = force->kspace->virial[1];
    kvirial_orig[2] = force->kspace->virial[2];
    kvirial_orig[3] = force->kspace->virial[3];
    kvirial_orig[4] = force->kspace->virial[4];
    kvirial_orig[5] = force->kspace->virial[5];
    
    if (update->eflag_atom) {
      double *keatom = force->kspace->eatom;
      for (i = 0; i < natom; i++)
        keatom_orig[i] = keatom[i];
    }
    if (update->vflag_atom) {
      double **kvatom = force->kspace->vatom;
      for (i = 0; i < natom; i++) {
        kvatom_orig[i][0] = kvatom[i][0];
        kvatom_orig[i][1] = kvatom[i][1];
        kvatom_orig[i][2] = kvatom[i][2];
        kvatom_orig[i][3] = kvatom[i][3];
        kvatom_orig[i][4] = kvatom[i][4];
        kvatom_orig[i][5] = kvatom[i][5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   restore arrays with force, energy, virial to original values
------------------------------------------------------------------------- */

void ComputeFEP::restore_accumulators()
{
  int i;

  int natom = atom->nlocal;

  double **f = atom->f;
  for (i = 0; i < natom; i++) {
    f[i][0] = f_orig[i][0]; 
    f[i][1] = f_orig[i][1]; 
    f[i][2] = f_orig[i][2]; 
  }

  force->pair->eng_vdwl = eng_vdwl_orig;
  force->pair->eng_coul = eng_coul_orig;

  force->pair->virial[0] = pvirial_orig[0];
  force->pair->virial[1] = pvirial_orig[1];
  force->pair->virial[2] = pvirial_orig[2];
  force->pair->virial[3] = pvirial_orig[3];
  force->pair->virial[4] = pvirial_orig[4];
  force->pair->virial[5] = pvirial_orig[5];

  if (update->eflag_atom) {
    double *peatom = force->pair->eatom;
    for (i = 0; i < natom; i++)
      peatom[i] = peatom_orig[i];
  }
  if (update->vflag_atom) {
    double **pvatom = force->pair->vatom;
    for (i = 0; i < natom; i++) {
      pvatom[i][0] = pvatom_orig[i][0];
      pvatom[i][1] = pvatom_orig[i][1];
      pvatom[i][2] = pvatom_orig[i][2];
      pvatom[i][3] = pvatom_orig[i][3];
      pvatom[i][4] = pvatom_orig[i][4];
      pvatom[i][5] = pvatom_orig[i][5];
    }
  }

  if (force->kspace) {
    force->kspace->energy = energy_orig;
    force->kspace->virial[0] = kvirial_orig[0];
    force->kspace->virial[1] = kvirial_orig[1];
    force->kspace->virial[2] = kvirial_orig[2];
    force->kspace->virial[3] = kvirial_orig[3];
    force->kspace->virial[4] = kvirial_orig[4];
    force->kspace->virial[5] = kvirial_orig[5];
    
    if (update->eflag_atom) {
      double *keatom = force->kspace->eatom;
      for (i = 0; i < natom; i++)
        keatom[i] = keatom_orig[i];
    }
    if (update->vflag_atom) {
      double **kvatom = force->kspace->vatom;
      for (i = 0; i < natom; i++) {
        kvatom[i][0] = kvatom_orig[i][0];
        kvatom[i][1] = kvatom_orig[i][1];
        kvatom[i][2] = kvatom_orig[i][2];
        kvatom[i][3] = kvatom_orig[i][3];
        kvatom[i][4] = kvatom_orig[i][4];
        kvatom[i][5] = kvatom_orig[i][5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   calculate pair energy of configuration using pair->single (test purposes)
------------------------------------------------------------------------- */

double ComputeFEP::energy_pair()
{
  int i,j,ii,jj,inum,jnum,itype,jtype,imol;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double fpair,factor_coul,factor_lj;
  double eng, eng_pair;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  pair = force->pair;
  double **cutsq = force->pair->cutsq;
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  fpair = 0.0;

  timer->stamp();
  eng = 0.0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;
      if (!(mask[j] & groupbit)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype])
        eng += pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
    }
  }
  timer->stamp(TIME_PAIR);

  eng_pair = 0.0;
  MPI_Allreduce(&eng,&eng_pair,1,MPI_DOUBLE,MPI_SUM,world);

  if (tailflag) {
    double volume = domain->xprd * domain->yprd * domain->zprd;
    eng_pair += force->pair->etail / volume; 
  }

  return eng_pair;
}
