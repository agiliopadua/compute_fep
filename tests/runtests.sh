#!/bin/bash

if [ $# == 0 ]; then
  echo usage: $0 command-to-run-lammps
  exit 1
fi

EXE=$@

. ./defs.sh

for tst in ${TSTSOFT} ${TSTADAPT} ${TSTFEP}; do
  echo ${tst}
  ${EXE} -in ${tst}.lmp -log ${tst}.log -var nsteps ${STEPS} \
     -echo none -screen none
done

