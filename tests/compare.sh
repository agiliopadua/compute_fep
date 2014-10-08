#!/bin/bash

. ./defs.sh

for tst in ${TSTSOFT} ${TSTADAPT}; do
  ${DIFFEXE} ${DIFFOPT} ${tst}.log ${REFDIR}/${tst}.log 
done

for tst in ${TSTFEP}; do
  ${DIFFEXE} ${DIFFOPT} ${tst}.log ${REFDIR}/${tst}.log 
  ${DIFFEXE} ${DIFFOPT} ${tst}.out ${REFDIR}/${tst}.out 
done
