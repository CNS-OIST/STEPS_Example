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
    partition, ind2sec = MorphPartition(mesh, morph, scale=1e-6)

    rois = []
    roi_areas = {}
    roi_vols = {}
    with mesh:
        for ind, tets in partition.tetTable.items():
            key = ind2sec[ind]
            roitet = ROI(tets, name=key)
            roitri = ROI(partition.triTable[ind], name=key + '_surf')
            rois.append(key)
            roi_vols[key] = roitet.Vol
            roi_areas[key] = roitri.Area

    return mesh, rois, roi_areas, roi_vols
