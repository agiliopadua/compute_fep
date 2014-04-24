#!/bin/bash

DIFFEXE=lmp_diff
DIFFOPT=numbers

TSTSOFT="lj lj00 ljsoft-lam00 ljsoft-lam05 ljsoft-lam10 ljsoft-ovrl-lam10"\
" tip4p tip4p00 tip4psoft-lam00 tip4psoft-lam05 tip4psoft-lam10"\
" tip4psoft-ovrl-lam10"

TSTFEP="fepNVT fepNPT"

for tst in ${TSTSOFT} ${TSTFEP}; do
  ${DIFFEXE} ${DIFFOPT} ${tst}.log ref/${tst}.log 
done

for tst in ${TSTFEP}; do
  ${DIFFEXE} ${DIFFOPT} ${tst}.log ref/${tst}.log 
  ${DIFFEXE} ${DIFFOPT} ${tst}.out ref/${tst}.out 
done
