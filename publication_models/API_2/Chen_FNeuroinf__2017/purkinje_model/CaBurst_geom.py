#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.interface

from steps.geom import *

########### MESH BRANCH MAPPING #################
def getGeom(mesh_file_name, morph_file_name = None):
    mesh = TetMesh.LoadAbaqus(mesh_file_name, scale=1e-06)

    with mesh:
        ########## Create an intracellular compartment i.e. cytosolic compartment
        cyto = TetComp.Create(mesh.tets, 'vsys')
        ########## Create a membrane as a surface mesh
        memb = TetPatch.Create(mesh.surface, cyto, None, 'ssys')
    
    if morph_file_name == None:
        return mesh
    
    morph = Morph.Load(morph_file_name)
    # partition based on branching
    partition = MorphPartition(mesh, morph, scale=1e-6)

    with mesh:
        for r, tets in partition.tetTable.items():
            ROI(tets, name=f'dend[{r}]')
        for r, tris in partition.triTable.items():
            ROI(tris, name=f'dend[{r}]_surf')

    return mesh
