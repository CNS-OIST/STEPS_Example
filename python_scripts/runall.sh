#!/bin/bash

# Usage:
#
# Run all tests: ./runall.sh
# Use another python executable: PYTHON=/path/to/another/python ./runall.sh
# Use another mpirun executable: MPIRUN=srun ./runall.sh
# Skip tests that require X: NO_X=1 ./run/all.sh

MPIRUN=${MPIRUN:-mpirun}
PYTHON=${PYTHON:-python}

${PYTHON} -c "import steps; steps._greet()"
echo "Example: Well Mixed"
(cd well_mixed && ${PYTHON} well_mixed.py)

echo "Example: IP3"
(cd ip3 && ${PYTHON} ip3.py)

echo "Example: Diffusion"
(cd diffusion && ${PYTHON} diffusion.py)

echo "Example: Surface Diffusion"
(cd surface_diffusion && ${PYTHON} surface_diffusion.py && ${PYTHON} surface_diffusion_tetode.py)

echo "Example: Diffusion Boundary"
(cd diffusion_boundary && ${PYTHON} diffusion_boundary.py)

echo "Example: Surface Diffusion Boundary"
(cd surface_diffusion_boundary && ${PYTHON} surface_diffusion_boundary.py && ${MPIRUN} -n 4 ${PYTHON} surface_diffusion_boundary_parallel.py)

echo "Example: SBML"
(cd sbml && ${PYTHON} sbml_import.py)

echo "Example: EField"
(cd HH_APprop && ${PYTHON} -c "import runHH_APprop; runHH_APprop.plotVz(10); runHH_APprop.plotVz(20); runHH_APprop.plotVz(30); runHH_APprop.show();  runHH_APprop.plotIz(10); runHH_APprop.plotIz(20); runHH_APprop.plotIz(30); runHH_APprop.show()" && ${MPIRUN} -n 4 ${PYTHON} HH_APprop_tetopsplit.py && ${PYTHON} -c "import runHH_APprop_tetode; runHH_APprop_tetode.plotVz(100); runHH_APprop_tetode.plotVz(200); runHH_APprop_tetode.plotVz(300); runHH_APprop_tetode.show()")

echo "Example: Parallelization, expect longer runtime"
(cd parallel && ${PYTHON} ssa_simulation.py && ${MPIRUN} -n 5 ${PYTHON} parallel_simulation.py && ${MPIRUN} -n 5 ${PYTHON} HH_APprop_tetopsplit.py)

echo "Example: Visualization Toolkit"
if [ "x${NO_X}" = x ] ;then
  cd visual && ${PYTHON} ip3r_sim.py)
fi
