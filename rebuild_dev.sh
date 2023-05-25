#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)

cd $REPO_ROOT
rm ${REPO_ROOT}/bunch1s.cpython-311-x86_64-linux-gnu.so || true
rm -rf ${REPO_ROOT}/build || true
mkdir build 
cmake -S . -B build 
cd build
make
cd $REPO_ROOT 
cp build/bunch1s.cpython-311-x86_64-linux-gnu.so ${REPO_ROOT}