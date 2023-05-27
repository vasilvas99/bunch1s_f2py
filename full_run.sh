#!/bin/bash
set -x

REPO_ROOT=$(git rev-parse --show-toplevel)
python3 ${REPO_ROOT}/integrate_trajectories.py
gfortran -O3 ${REPO_ROOT}/src/Fortran/MSI.f90 -o ${REPO_ROOT}/ms
${REPO_ROOT}/ms
python3 ${REPO_ROOT}/ms1_parser.py ${REPO_ROOT}/ms1_out/!!w-h__.dat
rm ${REPO_ROOT}/ms
