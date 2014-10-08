#!/bin/bash

STEPS=300

DIFFEXE=lmp_diff
DIFFOPT=numbers

REFDIR=ref-icms

TSTSOFT="lj lj00 ljsoft-lam00 ljsoft-lam05 ljsoft-lam10 ljsoft-ovrl-lam10"\
" tip4p tip4p00 tip4psoft-lam00 tip4psoft-lam05 tip4psoft-lam10"\
" tip4psoft-ovrl-lam10 charmm charmmsoft-lam10 ljsoft-lam10-omp"

TSTFEP="fepNVT fepNPT"
