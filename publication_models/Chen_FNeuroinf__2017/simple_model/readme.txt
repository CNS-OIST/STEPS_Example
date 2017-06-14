Scripts for parallel and SSA simulations of simple model in 

Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013

contact: w.chen@oist.jp

To record the time cost for each segments as described in the paper, please enable the add_definitions(-DMPI_PROFILING=1) flag in src/CMakeLists.txt in the source code before compiling STEPS.

General usage

Parallel simulation: 
mpirun -n NPROCS python parallel_simulation.py MOLECULE_RATIO MESHFILE RESULT_DIR

Serial SSA simulation:
python ssa_simulation.py MOLECULE_RATIO MESHFILE RESULT_DIR

NPROCS is the number of MPI processes used for simulation, NPROCS should be a multiple of 5, if not, adjust the following link in parallel_simulation.py:

tet_hosts = gd.linearPartition(g, [1,5,steps.mpi.nhosts/5])

to

tet_hosts = gd.linearPartition(g, [1,1,steps.mpi.nhosts])

Please note that this change will produce performance results different from the ones in the paper, due to the change of partitioning scheme.

MOLECULE_RATIO is ratio of molecule counts of current simulation against the default setting (simulations in Figure 4).
MESHFILE is the mesh file used in the simulation, all meshes are located in the “meshes” directory.
RESULT_DIR is the storage location for result data.

Specific usage for paper results:

Figure 3 (NPROCS starts from 50 to 300, with the increase of 50 each time):
mpirun -n NPROCS python parallel_simulation.py 1.0 10x10x100_3363tets.inp figure3


Figure 4 (NPROCS starts from 50 to 300, with the increase of 50 each time):
mpirun -n NPROCS python parallel_simulation.py 1.0 10x10x100_3363tets.inp figure4
mpirun -n NPROCS python parallel_simulation.py 0.1 10x10x100_3363tets.inp figure4
mpirun -n NPROCS python parallel_simulation.py 10.0 10x10x100_3363tets.inp figure4

python ssa_simulation.py 1.0 10x10x100_3363tets.inp figure4
python ssa_simulation.py 0.1 10x10x100_3363tets.inp figure4
python ssa_simulation.py 10.0 10x10x100_3363tets.inp figure4

Figure 5 (NPROCS starts from 50 to 300, with the increase of 50 each time): 
mpirun -n NPROCS python parallel_simulation.py 1.0 10x10x100_3363tets.inp figure5
mpirun -n NPROCS python parallel_simulation.py 1.0 10x10x100_13009tets.inp figure5
mpirun -n NPROCS python parallel_simulation.py 1.0 10x10x100_113096tets.inp figure5

Figure 6:
mpirun -n 300 python parallel_simulation.py 1.0 10x10x100_3363tets.inp figure6
mpirun -n 300 python parallel_simulation.py 1.0 10x10x200.inp figure6
mpirun -n 300 python parallel_simulation.py 1.0 10x10x300.inp figure6
mpirun -n 300 python parallel_simulation.py 1.0 10x20x100.inp figure6
mpirun -n 300 python parallel_simulation.py 1.0 10x30x100.inp figure6