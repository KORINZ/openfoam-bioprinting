#!/bin/sh
#rm -r constant/polyMesh/*
#cp -r /home/colin/OpenFOAM/colin-9/run/wsl/mesh/current_needle_mesh/* constant/polyMesh
renumberMesh -overwrite
#transformPoints -scale '(0.001 0.001 0.001)'
decomposePar -force
mpirun -np 6 pimpleFoam -parallel 1> foam.log 2>> error.log &&
reconstructPar > ReconstructPar.log &&
rm -r processor* &&
foamLog foam.log &
#simpleFoam -postProcess -func shearStress &&
