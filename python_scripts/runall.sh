#!/bin/bash
python -c "import steps; steps._greet()"
echo "Example: Well Mixed"
cd well_mixed && python well_mixed.py && cd ..
echo "Example: IP3"
cd ip3 && python ip3.py && cd ..
echo "Example: Diffusion"
cd diffusion && python diffusion.py && cd ..
echo "Example: Surface Diffusion"
cd surface_diffusion && python surface_diffusion.py && python surface_diffusion_tetode.py && cd ..
echo "Example: Diffusion Boundary"
cd diffusion_boundary && python diffusion_boundary.py && cd ..
echo "Example: Surface Diffusion Boundary"
cd surface_diffusion_boundary && python surface_diffusion_boundary.py && mpirun -n 4 python surface_diffusion_boundary_parallel.py && cd ..
echo "Example: SBML"
cd sbml && python sbml_import.py && cd ..
echo "Example: EField"
cd HH_APprop && python -c "import runHH_APprop; runHH_APprop.plotVz(10); runHH_APprop.plotVz(20); runHH_APprop.plotVz(30); runHH_APprop.show();  runHH_APprop.plotIz(10); runHH_APprop.plotIz(20); runHH_APprop.plotIz(30); runHH_APprop.show()" && mpirun -n 4 python HH_APprop_tetopsplit.py && python -c "import runHH_APprop_tetode; runHH_APprop_tetode.plotVz(100); runHH_APprop_tetode.plotVz(200); runHH_APprop_tetode.plotVz(300); runHH_APprop_tetode.show()" && cd ..
echo "Example: Parallelization, expect longer runtime"
cd parallel && python ssa_simulation.py && mpirun -n 5 python parallel_simulation.py && mpirun -n 5 python HH_APprop_tetopsplit.py && cd ..
echo "Example: Visualization Toolkit"
cd visual && python ip3r_sim.py && cd ..