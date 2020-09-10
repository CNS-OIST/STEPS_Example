#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

from steps.geom import *

import os
from subprocess import call

os.makedirs("meshes/partition", exist_ok=True)

MESH_FILE = "meshes/branch.inp"
mesh = TetMesh.LoadAbaqus(MESH_FILE, scale=1e-06)

mesh.ConvertToMetis('meshes/partition/branch.metis')

print("Generate partition for desktop computer (from 2 cores to 10 cores)")
for i in range(2, 11, 2):
    call(['mpmetis', '-ncommon=3', '-minconn', '-niter=1000', 'meshes/partition/branch.metis', '%i' % (i)])
    metis_part = MetisPartition(mesh, 'meshes/partition/branch.metis.epart.%i' % (i), default_tris=mesh.surface)
    metis_part.printStats()

print("Generate partition for supercomputer (from 50 cores to 1000 cores)")
for i in range(50, 1001, 50):
    call(['mpmetis', '-ncommon=3', '-minconn', '-niter=1000', 'meshes/partition/branch.metis', '%i' % (i)])
    metis_part = MetisPartition(mesh, 'meshes/partition/branch.metis.epart.%i' % (i), default_tris=mesh.surface)
    metis_part.printStats()
