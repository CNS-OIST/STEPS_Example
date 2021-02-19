Scripts for parallel and SSA simulations of Purkinje cell sub-branch model with calcium influx in 

Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013

contact: w.chen@oist.jp

To record the time cost for each segments as described in the paper, please enable the add_definitions(-DMPI_PROFILING=1) flag in src/CMakeLists.txt in the source code before compiling STEPS.

Specific usage for paper results

Figure 7 (Rquire PyQtGraph and PyOpenGL installed):

python view_partitions.py meshes/partition/branch.metis.epart.50
python view_partitions.py meshes/partition/branch.metis.epart.1000

Figure 8 (NPROCS starts from 50 to 1000, with the increase of 50 each time):

mpirun -n NPROCS python parallel_branch_influx.py figure8_influx
mpirun -n NPROCS python parallel_branch_background.py figure8_background
python ssa_branch_influx.py figure8_influx
python ssa_branch_influx.py figure8_background

If you want to visualize the recorded calcium activity such as the one shown in Figure 8A, you need to first install PyQtGraph, PyOpenGL, NumPy and Matplotlib, and use the following command:

python view_activity.py YOUR_SIM_RECORDING

An example recording file is provided in the ca_activity_result_example directory.

Figure 9 (require add_definitions(-DMPI_PROFILING=1) flag in src/CMakeLists.txt enable):

mpirun -n 50 python parallel_branch_influx.py figure9

Figure 10 (NPROCS starts from 100 to 2000, with the increase of 100 each time):

mpirun -n NPROCS python parallel_fullcell_background.py figure10

Figure 11:

Simulations in Figure 8, 10, and the following 

python ssa_fullcell_background.py figure11

Other General Usage:

1. To convert a NEUERON hoc cell morphology to STEPPS morph structure (Required in influx simulation for branch mapping):
python neuron2morph.py meshes/branch.hoc meshes/branch.morph

2. To partition the branch or full cell morphology using Metis:
python branch_partition.py
python fullcell_partition.py

partitions can be visualized using the visualization toolkit:

python view_partitions.py MESH_FILE PARTITION_FILE, e.g.

python view_partitions.py meshes/branch.inp meshes/partition/branch.metis.epart.1000
python view_partitions.py meshes/fullcell.inp meshes/partition/fullcell.metis.epart.1000


