.. _API_2_geom:

****************
steps.API_2.geom
****************

The ``geom`` module contains classes related to meshes and geometries: compartments, patches, region of interest, diffusion boundaries, etc.
It also contains utility classes for handling element lists.

======================
Detailed documentation
======================

* Containers

 * :py:class:`steps.API_2.geom.Geometry`
 * :py:class:`steps.API_2.geom.TetMesh`
 * :py:class:`steps.API_2.geom.DistMesh`

* Geometrical elements

  * :py:class:`steps.API_2.geom.Compartment`
  * :py:class:`steps.API_2.geom.Patch`
  * :py:class:`steps.API_2.geom.Membrane`
  * :py:class:`steps.API_2.geom.ROI`
  * :py:class:`steps.API_2.geom.DiffBoundary`
  * :py:class:`steps.API_2.geom.SDiffBoundary`

* Mesh partitioning

 * :py:class:`steps.API_2.geom.MeshPartition`
 * :py:class:`steps.API_2.geom.LinearMeshPartition`
 * :py:class:`steps.API_2.geom.MetisPartition`
 * :py:class:`steps.API_2.geom.MorphPartition`

* References

 * :py:class:`steps.API_2.geom.Reference`
 * :py:class:`steps.API_2.geom.TetReference`
 * :py:class:`steps.API_2.geom.TriReference`
 * :py:class:`steps.API_2.geom.BarReference`
 * :py:class:`steps.API_2.geom.VertReference`

* Reference lists

 * :py:class:`steps.API_2.geom.RefList`
 * :py:class:`steps.API_2.geom.TetList`
 * :py:class:`steps.API_2.geom.TriList`
 * :py:class:`steps.API_2.geom.BarList`
 * :py:class:`steps.API_2.geom.VertList`

* Convenience classes

 * :py:class:`steps.API_2.geom.Point`
 * :py:class:`steps.API_2.geom.BoundingBox`
 * :py:class:`steps.API_2.geom.Morph`

----------

.. automodule:: steps.API_2.geom
    :members:
    :undoc-members:
    :inherited-members: ndarray
    :special-members:

