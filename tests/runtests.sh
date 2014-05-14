#!/bin/bash

if [ $# == 0 ]; then
  echo usage: $0 command-to-run-lammps
  exit 1
fi

EXE=$@
STEPS=300

TSTSOFT="lj lj00 ljsoft-lam00 ljsoft-lam05 ljsoft-lam10 ljsoft-ovrl-lam10"\
" tip4p tip4p00 tip4psoft-lam00 tip4psoft-lam05 tip4psoft-lam10"\
" tip4psoft-ovrl-lam10"

TSTFEP="fepNVT fepNPT"

for tst in ${TSTSOFT}; do
  echo ${tst}
  ${EXE} -in ${tst}.in -log ${tst}.log -var nsteps ${STEPS} \
     -echo none -screen none
done

for tst in ${TSTFEP}; do
  echo ${tst}
  ${EXE} -in ${tst}.in -log ${tst}.log -var nsteps ${STEPS} \
     -echo none -screen none
done
