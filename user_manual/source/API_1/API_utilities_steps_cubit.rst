.. _API_1_utilities_steps_cubit:

*********************************************************************
API Reference: STEPS-CUBIT Geometry Preparation Toolkit (Independent)
*********************************************************************

Shadow Classes
==============
* :class:`steps.API_1.utilities.steps_cubit.ShadowMesh`
* :class:`steps.API_1.utilities.steps_cubit.ShadowComp`
* :class:`steps.API_1.utilities.steps_cubit.ShadowPatch`

Functions
=========
.. automodule:: steps.API_1.utilities.steps_cubit
   :members: ShadowMesh, ShadowComp, ShadowPatch, getSelectedVolumes, getSelectedSurfaces, getSelectedNodes, getSelectedTets, getSelectedTris, selectedVolumesAsComp, selectedSurfacesAsPatch, selectedTetsAsComp, selectedTrisAsPatch, selectedNodesAsROI, selectedTetsAsROI, selectedTrisAsROI, getNodesBoundInSelectedVols, getTetsBoundInSelectedVols, getTrisBoundInSelectedVols, boundTetsAsComp, boundTrisAsPatch, boundNodesAsROI, boundTetsAsROI, boundTrisAsROI, drawROI, highlightROI, toStr
.. autofunction:: drawComp(mesh, comp_id)
.. autofunction:: drawPatch(mesh, patch_id)
.. autofunction:: highlightComp(mesh, comp_id)
.. autofunction:: highlightPatch(mesh, patch_id)
