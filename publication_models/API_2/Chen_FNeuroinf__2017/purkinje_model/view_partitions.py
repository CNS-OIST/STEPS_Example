#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

from steps.geom import *
from steps.visual import *

import sys

try:
    _, MESH_FILE, PARTITION_FILE  = sys.argv
except:
    MESH_FILE = "meshes/branch.inp"
    PARTITION_FILE = "meshes/partition/branch.metis.epart.100"

mesh = TetMesh.LoadAbaqus(MESH_FILE, scale=1e-06)
partition = MetisPartition(mesh, PARTITION_FILE)

sc = SimControl()
with sc:
    PartitionDisplay(partition)
sc.run()

