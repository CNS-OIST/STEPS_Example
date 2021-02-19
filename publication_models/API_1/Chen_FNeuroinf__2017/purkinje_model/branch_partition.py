#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

from __future__ import print_function
import steps.utilities.geom_decompose as gd
import steps.utilities.meshio as meshio
import steps.utilities.metis_support as metis
import os

try: os.mkdir("meshes/partition")
except: pass

MESH_FILE = "meshes/branch.inp"
mesh = meshio.importAbaqus(MESH_FILE, 1e-6)[0]

metis.tetmesh2metis(mesh, 'meshes/partition/branch.metis')

from subprocess import call
print("Generate partition for desktop computer (from 2 cores to 10 cores)")
for i in range(2, 11, 2):
    call(['mpmetis', '-ncommon=3', '-minconn', '-niter=1000', 'meshes/partition/branch.metis', '%i' % (i)])
    metis_parts = metis.readPartition('meshes/partition/branch.metis.epart.%i' % (i))
    surf_tris = mesh.getSurfTris()
    tri_parts = gd.partitionTris(mesh, metis_parts, surf_tris)
    gd.printPartitionStat(metis_parts, tri_parts)

print("Generate partition for supercomputer (from 50 cores to 1000 cores)")
for i in range(50, 1001, 50):
    call(['mpmetis', '-ncommon=3', '-minconn', '-niter=1000', 'meshes/partition/branch.metis', '%i' % (i)])
    metis_parts = metis.readPartition('meshes/partition/branch.metis.epart.%i' % (i))
    surf_tris = mesh.getSurfTris()
    tri_parts = gd.partitionTris(mesh, metis_parts, surf_tris)
    gd.printPartitionStat(metis_parts, tri_parts)
