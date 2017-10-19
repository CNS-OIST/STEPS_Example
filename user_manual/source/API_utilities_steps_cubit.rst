.. _API_utilities_steps_cubit:

*********************************************************************
API Reference: STEPS-CUBIT Geometry Preparation Toolkit (Independent)
*********************************************************************

Shadow Classes
==============
* :class:`steps_cubit.ShadowMesh`
* :class:`steps_cubit.ShadowComp`
* :class:`steps_cubit.ShadowPatch`

.. automodule:: steps_cubit
   :members: ShadowMesh, ShadowComp, ShadowPatch

Functions
=========
.. automodule:: steps_cubit
   :members: getSelectedVolumes, getSelectedSurfaces, getSelectedNodes, getSelectedTets, getSelectedTris, selectedVolumesAsComp, selectedSurfacesAsPatch, selectedTetsAsComp, selectedTrisAsPatch, selectedNodesAsROI, selectedTetsAsROI, selectedTrisAsROI, getNodesBoundInSelectedVols, getTetsBoundInSelectedVols, getTrisBoundInSelectedVols, boundTetsAsComp, boundTrisAsPatch, boundNodesAsROI, boundTetsAsROI, boundTrisAsROI, drawROI, highlightROI, toStr
.. autofunction:: drawComp(comp)
.. autofunction:: drawComp(mesh, comp_id)
.. autofunction:: drawPatch(patch)
.. autofunction:: drawPatch(mesh, patch_id)
.. autofunction:: highlightComp(comp)
.. autofunction:: highlightComp(mesh, comp_id)
.. autofunction:: highlightPatch(patch)
.. autofunction:: highlightPatch(mesh, patch_id)
