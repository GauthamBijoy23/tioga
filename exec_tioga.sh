#!/usr/bin/env bash
cd ~/tioga/build
make
cd ~/tioga/case
mpirun -np 1 ~/tioga/build/driver/tioga_read.exe
cd ~/tioga/build/driver/output/


