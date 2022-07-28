.. _interfaces:

**************************
Python interfaces to STEPS
**************************

Since version 3.6, STEPS is shipped with a new, more pythonic interface. In this documentation, we will refer to the original interface as API_1 and to the new one as API_2. Unless mentioned explicitely, we will use API_2 by default.

Compatibility between APIs
==========================

Models that were written using API_1 are still usable without any changes. Users are however encouraged to start using API_2 for new models. It is not advised to try to mix both APIs in a single python script.
Until API_1 becomes deprecated, users that want to use API_2 have to explicitely specify it at the beginning of their python scripts (see :doc:`API_2/API_ref`).

Guides
======

Chapters 3 to 9 introduce STEPS concepts in details through examples:

    - :doc:`API_2/Interface_Tutorial_1_wm`
    - :doc:`API_2/Interface_Tutorial_2_IP3`
    - :doc:`API_2/Interface_Tutorial_3_Diffusion`
    - :doc:`API_2/Interface_Tutorial_4_Complexes`
    - :doc:`API_2/Interface_Tutorial_5_Efield`
    - :doc:`API_2/Interface_Tutorial_6_MPI`
    - :doc:`API_2/Interface_Tutorial_7_visual`
    - :doc:`API_2/Interface_Tutorial_8_CaBurst`

Equivalent chapters using API_1 are still available for reference (see :doc:`API_1_guide`).

Finally, :doc:`API_2/API_ref` provides detailed documentation of STEPS API, documentation relative to API_1 is still available (:doc:`API_1/API_ref`).

All chapters from the guide are runnable with [jupyter lab](https://jupyter.org/). The corresponding `*.ipynb` files are
available in the [STEPS_Example](https://github.com/CNS-OIST/STEPS_Example/tree/master/user_manual/source/API_2) github repository.

In the chapters, jupyter code cells are represented with e.g.:

.. raw:: html

    <div class="nbinput nblast docutils container">
    <div class="prompt highlight-none notranslate"><div class="highlight"><pre><span></span>[1]:
    </pre></div>
    </div>
    <div class="input_area highlight-ipython3 notranslate"><div class="highlight"><pre><span></span><span class="kn">import</span> <span class="nn">steps.interface</span>

    <span class="kn">from</span> <span class="nn">steps.model</span> <span class="kn">import</span> <span class="o">*</span>
    <span class="kn">from</span> <span class="nn">steps.geom</span> <span class="kn">import</span> <span class="o">*</span>
    <span class="kn">from</span> <span class="nn">steps.rng</span> <span class="kn">import</span> <span class="o">*</span>
    <span class="kn">from</span> <span class="nn">steps.sim</span> <span class="kn">import</span> <span class="o">*</span>
    <span class="kn">from</span> <span class="nn">steps.saving</span> <span class="kn">import</span> <span class="o">*</span>

    <span class="kn">from</span> <span class="nn">matplotlib</span> <span class="kn">import</span> <span class="n">pyplot</span> <span class="k">as</span> <span class="n">plt</span>
    <span class="kn">import</span> <span class="nn">numpy</span> <span class="k">as</span> <span class="nn">np</span>

    <span class="n">mdl</span> <span class="o">=</span> <span class="n">Model</span><span class="p">()</span>

    <span class="n">r</span> <span class="o">=</span> <span class="n">ReactionManager</span><span class="p">()</span>
    </pre></div>
    </div>
    </div>
