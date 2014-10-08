#!/bin/bash

if [ $# == 0 ]; then
  echo usage: $0 command-to-run-lammps
  exit 1
fi

EXE=$@

. ./defs.sh

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
