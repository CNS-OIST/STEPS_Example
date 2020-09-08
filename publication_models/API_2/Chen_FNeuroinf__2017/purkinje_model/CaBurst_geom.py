import steps.interface

#########################################################################
#  This script is provided for
#
#  Chen W and De Schutter E (2017) Parallel STEPS: Large Scale Stochastic Spatial Reaction-Diffusion Simulation with High Performance Computers. Front. Neuroinform. 11:13. doi: 10.3389/fninf.2017.00013
#
##########################################################################

from steps.geom import *
import steps.utilities.morph_support as morph_support
try:
    import cPickle as pickle
except:
    import pickle

########### MESH BRANCH MAPPING #################
def getGeom(mesh_file_name, morph_file_name = None):
    mesh = TetMesh.LoadAbaqus(mesh_file_name, scale=1e-06)

    ########## Create an intracellular compartment i.e. cytosolic compartment
    inner_tets = range(len(mesh.tets))
    with mesh:
        cyto = TetComp.Create(inner_tets, 'vsys')

        ########## Create a membrane as a surface mesh
        surf_tris = mesh.surface.indices
        memb = TetPatch.Create(surf_tris, cyto, None, 'ssys')
    
    if morph_file_name == None:
        return mesh
    
    # morph file is a cPickled dictionary of branching data from NEURON .hoc file, neuron2morph.py for detail
    morph_file = open(morph_file_name, 'r')
    morph = pickle.load(morph_file)
    
    # partition tetrahedrons based on branching
    branch_tets = morph_support.mapMorphTetmesh(morph, mesh)

    # partition surface triangles based on above tetrahedron partition
    # WARNING: partitionTris was incorporated into LinearMeshPartition or MetisPartition.
    ...

    branch_tet_table = gd.getTetPartitionTable(branch_tets)
    branch_tri_table = gd.getTriPartitionTable(branch_surf_tris)

    rois = []
    roi_areas = {}
    roi_vols = {}
    
    # add the branch mapping as ROI
    for r in range(101):
        # WARNING: Using a variable name that is reserved (['r']).
        roi = 'dend[%i]' % (r)
        with mesh:
            ROI(TetList(branch_tet_table[roi]), name=roi)
        with mesh:
            ROI(TriList(branch_tri_table[roi]), name="%s_surf" % (roi))
        rois.append(roi)
        roi_areas[roi] = getattr(mesh, "%s_surf" % (roi)).Area
        roi_vols[roi] = getattr(mesh, roi).Vol
    
    return mesh, rois, roi_areas, roi_vols
