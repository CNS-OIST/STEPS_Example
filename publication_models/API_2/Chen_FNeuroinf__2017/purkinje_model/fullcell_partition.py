import steps.interface

#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

from __future__ import print_function
import os

try: os.mkdir("meshes/partition")
except: pass

MESH_FILE = "meshes/fullcell.inp"
mesh = TetMesh.LoadAbaqus(MESH_FILE, scale=1e-06)

metis.tetmesh2metis(mesh, 'meshes/partition/fullcell.metis')

from subprocess import call
print("Generate partition for desktop computer (from 2 cores to 10 cores)")
for i in range(2, 11, 2):
    call(['mpmetis', '-ncommon=3', '-minconn', '-niter=1000', 'meshes/partition/fullcell.metis', '%i' % (i)])
    metis_parts = metis.readPartition('meshes/partition/fullcell.metis.epart.%i' % (i))
    surf_tris = mesh.surface.indices
    # WARNING: partitionTris was incorporated into LinearMeshPartition or MetisPartition.
    ...
    gd.printPartitionStat(metis_parts, tri_parts)

print("Generate partition for supercomputer (from 100 cores to 2000 cores)")
for i in range(2100, 5001, 100):
    call(['mpmetis', '-ncommon=3', '-minconn', '-niter=1000', 'meshes/partition/fullcell.metis', '%i' % (i)])
    metis_parts = metis.readPartition('meshes/partition/fullcell.metis.epart.%i' % (i))
    surf_tris = mesh.surface.indices
    # WARNING: partitionTris was incorporated into LinearMeshPartition or MetisPartition.
    ...
    gd.printPartitionStat(metis_parts, tri_parts)
