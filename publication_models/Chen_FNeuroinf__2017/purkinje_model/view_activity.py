#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

from extra.activity_viewer import *
import cPickle
import steps.utilities.geom_decompose as gd
import steps.utilities.meshio as meshio
import steps.visual
import pyqtgraph as pg
import random
import math
import sys

ACTIVITY_FILE = sys.argv[1]

MESH_FILE = "meshes/branch.inp"
mesh = meshio.importAbaqus(MESH_FILE, 1e-6)[0]
morph_file = open("meshes/branch.morph", 'r')
morph = cPickle.load(morph_file)
tet_parts = gd.mapMorphTetmesh(morph, mesh)
tri_parts = gd.partitionTris(mesh, tet_parts, mesh.getSurfTris())

tet_part_table = gd.getTetPartitionTable(tet_parts)
tri_part_table = gd.getTriPartitionTable(tri_parts)

mpi_data = SI2NEURON(readData(ACTIVITY_FILE))
app = pg.mkQApp()
conc_display = TriActivitySeriesDisplay(mesh, tri_parts, mpi_data, title = "MPI Sim", time_unit = "ms", data_unit = "mM", min_v = 0.0236, max_v = 3.51)
app.exec_()
