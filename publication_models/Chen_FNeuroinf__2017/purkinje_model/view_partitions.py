#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

try:
    import cPickle as pickle
except:
    import pickle
import steps.utilities.geom_decompose as gd
import steps.utilities.meshio as meshio
import steps.utilities.metis_support as metis
import sys

if len(sys.argv) == 3:
    MESH_FILE = sys.argv[1]
    PARTITION_FILE = sys.argv[2]
else:
    MESH_FILE = "meshes/branch.inp"
    PARTITION_FILE = "meshes/partition/branch.metis.epart.100"

mesh = meshio.importAbaqus(MESH_FILE, 1e-6)[0]
tet_partitions = metis.readPartition(PARTITION_FILE)
tri_partitions = gd.partitionTris(mesh, tet_partitions, mesh.getSurfTris())

import steps.visual
import pyqtgraph as pg
app = pg.mkQApp()

w = steps.visual.TriPartitionDisplay(mesh, tri_partitions, w = 1200, h = 800)
app.exec_()
