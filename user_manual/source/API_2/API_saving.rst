.. _API_2_saving:

******************
steps.API_2.saving
******************

The ``saving`` module contains classes related to data saving during simulation.

======================
Detailed documentation
======================

* Result selectors

    * :py:class:`steps.API_2.saving.ResultSelector`
    * :py:class:`steps.API_2.saving.CustomResults`

* Saving to database

    * :py:class:`steps.API_2.saving.DatabaseHandler`

        * :py:class:`steps.API_2.saving.SQLiteDBHandler`
        * :py:class:`steps.API_2.saving.HDF5Handler`
        * :py:class:`steps.API_2.saving.XDMFHandler`

    * :py:class:`steps.API_2.saving.DatabaseGroup`

        * :py:class:`steps.API_2.saving.HDF5Group`
        * :py:class:`steps.API_2.saving.SQLiteGroup`

* Reading from database

    * :py:class:`steps.API_2.saving.HDF5Handler`
    * :py:class:`steps.API_2.saving.HDF5MultiFileReader`

----------

.. automodule:: steps.API_2.saving
    :members:
    :undoc-members:
    :special-members:
    :show-inheritance:


