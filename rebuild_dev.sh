#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)

cd $REPO_ROOT
find . -name "bunch1s.*.so" -delete
rm -rf ${REPO_ROOT}/build || true
mkdir build 
cmake -S . -B build 
cd build
make
cd $REPO_ROOT 
cp build/bunch1s.*.so ${REPO_ROOT}