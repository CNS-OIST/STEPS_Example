.. _API_2_model:

*****************
steps.API_2.model
*****************

The ``model`` module contains classes related to model definition: species, complexes, reactions, diffusion rules, channels, etc.

======================
Detailed documentation
======================

* Containers

 * :py:class:`steps.API_2.model.Model`
 * :py:class:`steps.API_2.model.VolumeSystem`
 * :py:class:`steps.API_2.model.SurfaceSystem`

* Reactants

 * :py:class:`steps.API_2.model.Species`
 * :py:class:`steps.API_2.model.Complex`

  * :py:class:`steps.API_2.model.ComplexState`
  * :py:class:`steps.API_2.model.ComplexSelector`
  * :py:class:`steps.API_2.model.SubUnit`

   * :py:class:`steps.API_2.model.SubUnitState`
   * :py:class:`steps.API_2.model.SubUnitSelector`
 
  * :py:func:`steps.API_2.model.NoOrdering`
  * :py:func:`steps.API_2.model.StrongOrdering`
  * :py:func:`steps.API_2.model.RotationalSymmetryOrdering`
 
 * :py:class:`steps.API_2.model.Channel`

* Reaction-diffusion

 * :py:class:`steps.API_2.model.ReactionManager`
 * :py:class:`steps.API_2.model.Reaction`
 
  * :py:class:`steps.API_2.model.VDepRate`
  * :py:class:`steps.API_2.model.CompDepRate`
 
 * :py:class:`steps.API_2.model.Diffusion`

  * :py:class:`steps.API_2.model.CompDepDcst`

* Currents

 * :py:class:`steps.API_2.model.OhmicCurr`

  * :py:class:`steps.API_2.model.CompDepCond`
 
 * :py:class:`steps.API_2.model.GHKCurr`

  * :py:class:`steps.API_2.model.CompDepP`
  * :py:func:`steps.API_2.model.GHKCurr.PInfo`

----------

.. automodule:: steps.API_2.model
    :members:
    :undoc-members:
    :inherited-members: ndarray
    :special-members:

