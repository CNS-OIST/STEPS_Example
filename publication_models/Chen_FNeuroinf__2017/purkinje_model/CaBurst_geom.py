#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

import steps.geom as sgeom
import steps.utilities.meshio as meshio
import steps.utilities.geom_decompose as gd
import steps.utilities.morph_support as morph_support
try:
    import cPickle as pickle
except:
    import pickle

########### MESH BRANCH MAPPING #################
def getGeom(mesh_file_name, morph_file_name = None):
    mesh = meshio.importAbaqus(mesh_file_name, 1e-6)[0]

    ########## Create an intracellular compartment i.e. cytosolic compartment
    inner_tets = range(mesh.ntets)
    cyto = sgeom.TmComp('cyto', mesh, inner_tets)
    cyto.addVolsys('vsys')

    ########## Create a membrane as a surface mesh
    surf_tris = mesh.getSurfTris()
    memb = sgeom.TmPatch('memb', mesh, surf_tris, cyto)
    memb.addSurfsys('ssys')
    
    if morph_file_name == None: return mesh
    
    # morph file is a cPickled dictionary of branching data from NEURON .hoc file, neuron2morph.py for detail
    morph_file = open(morph_file_name, 'r')
    morph = pickle.load(morph_file)
    
    # partition tetrahedrons based on branching
    branch_tets = morph_support.mapMorphTetmesh(morph, mesh)

    # partition surface triangles based on above tetrahedron partition
    branch_surf_tris = gd.partitionTris(mesh, branch_tets, surf_tris)

    branch_tet_table = gd.getTetPartitionTable(branch_tets)
    branch_tri_table = gd.getTriPartitionTable(branch_surf_tris)

    rois = []
    roi_areas = {}
    roi_vols = {}
    
    # add the branch mapping as ROI
    for r in range(101):
        roi = 'dend[%i]' % (r)
        mesh.addROI(roi, sgeom.ELEM_TET, branch_tet_table[roi])
        mesh.addROI("%s_surf" % (roi), sgeom.ELEM_TRI, branch_tri_table[roi])
        rois.append(roi)
        roi_areas[roi] = mesh.getROIArea("%s_surf" % (roi))
        roi_vols[roi] = mesh.getROIVol(roi)
    
    return mesh, rois, roi_areas, roi_vols
