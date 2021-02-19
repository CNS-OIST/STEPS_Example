.. _API_1_mpisolver:

**********************
steps.API_1.mpi.solver
**********************

.. data:: steps.API_1.mpi.nhosts

The number of processes used in the MPI simulation.

.. data:: steps.API_1.mpi.rank

The rank of the current MPI process.

.. module:: steps.API_1.mpi.solver

Implementation of parallel simulation solver.
   
* :class:`steps.API_1.mpi.solver.TetOpSplit`

The solver is a partial implementation of the STEPS solver API.
At the moment STEPS implements just one parallel solver. 

:class:`steps.API_1.mpi.solver.TetOpSplit` implements a spatial stochastic 
solver that applies an accurate approximation for diffusion, and reactions
are simulated stochastically based on Gillespie's Direct SSA Method. Simulations
are run in parallel under the MPI protocol.

Several EField solution are provided in the current parallel implementation.
For large scale parallel simulation with hundreds of processes we recommand `EF_DV_PETSC`.
For small scale simulation with several process we recommand `EF_DV_BDSYS`.

.. autoclass:: TetOpSplit

    Options for `calcMembPot`

    * .. data:: steps.API_1.mpi.solver.EF_NONE
    
         No EField solver is needed.

    * .. data:: steps.API_1.mpi.solver.EF_DEFAULT

         Run serial EField simulation (Tetexact version) on process 0.

    * .. data:: steps.API_1.mpi.solver.EF_DV_BDSYS

         Use parallel SuperLU EField solver.

    * .. data:: steps.API_1.mpi.solver.EF_DV_PETSC

         Use parallel PETSc EField solver.

    **Solver Information**
    
    .. automethod:: getSolverName
    .. automethod:: getSolverDesc
    .. automethod:: getSolverAuthors
    .. automethod:: getSolverEmail
    
    **Solver Control**
    
    .. automethod:: reset
    .. automethod:: checkpoint
    .. automethod:: restore
    .. automethod:: run
    .. automethod:: advance
    .. automethod:: step
    .. automethod:: getTime
    .. automethod:: getA0
    .. automethod:: getNSteps
    .. automethod:: getTemp
    .. automethod:: setTemp
    .. automethod:: setEfieldDT
    .. automethod:: getCompTime
    .. automethod:: getSyncTime
    .. automethod:: getIdleTime
    .. automethod:: getReacExtent
    .. automethod:: getDiffExtent
    .. automethod:: getNIteration
    .. automethod:: repartitionAndReset
    
    **Compartment Data Access**
    
    .. automethod:: getCompVol
    .. automethod:: getCompCount
    .. automethod:: setCompCount
    .. automethod:: getCompAmount
    .. automethod:: setCompAmount
    .. automethod:: getCompConc
    .. automethod:: setCompConc
    .. automethod:: getCompClamped
    .. automethod:: setCompClamped 
    .. automethod:: getCompReacK 
    .. automethod:: setCompReacK
    .. automethod:: getCompReacActive
    .. automethod:: setCompReacActive
    .. automethod:: getCompDiffD
    .. automethod:: setCompDiffD
    .. automethod:: getCompDiffActive
    .. automethod:: setCompDiffActive
    .. automethod:: getCompReacC
    .. automethod:: getCompReacH
    .. automethod:: getCompReacA
    .. automethod:: getCompReacExtent
    .. automethod:: resetCompReacExtent
    
    **Patch Data Access**
    
    .. automethod:: getPatchArea
    .. automethod:: getPatchCount
    .. automethod:: setPatchCount
    .. automethod:: getPatchAmount
    .. automethod:: setPatchAmount
    .. automethod:: getPatchClamped
    .. automethod:: setPatchClamped
    .. automethod:: getPatchSReacK
    .. automethod:: setPatchSReacK
    .. automethod:: getPatchSReacActive
    .. automethod:: setPatchSReacActive
    .. automethod:: getPatchSReacC
    .. automethod:: getPatchSReacH
    .. automethod:: getPatchSReacA
    .. automethod:: getPatchSReacExtent
    .. automethod:: resetPatchSReacExtent
    .. automethod:: getPatchVDepSReacActive
    .. automethod:: setPatchVDepSReacActive
    
    **Diffusion Boundary Data Access**
    
    .. automethod:: setDiffBoundaryDiffusionActive
    .. automethod:: getDiffBoundaryDiffusionActive
    .. automethod:: setDiffBoundaryDcst
    
    **Surface Diffusion Boundary Data Access**
    
    .. automethod:: setSDiffBoundaryDiffusionActive
    .. automethod:: getSDiffBoundaryDiffusionActive
    .. automethod:: setSDiffBoundaryDcst

    **Tetrahedral Data Access**
    
    .. automethod:: getTetVol
    .. automethod:: getTetSpecDefined
    .. automethod:: getTetCount
    .. automethod:: setTetCount
    .. automethod:: getTetAmount
    .. automethod:: setTetAmount
    .. automethod:: getTetConc
    .. automethod:: setTetConc
    .. automethod:: getTetClamped
    .. automethod:: setTetClamped
    .. automethod:: getTetReacK
    .. automethod:: setTetReacK
    .. automethod:: getTetReacActive
    .. automethod:: setTetReacActive
    .. automethod:: getTetDiffD
    .. automethod:: setTetDiffD
    .. automethod:: getTetDiffActive
    .. automethod:: setTetDiffActive
    .. automethod:: getTetReacC
    .. automethod:: getTetReacH
    .. automethod:: getTetReacA
    .. automethod:: getTetDiffA
    .. automethod:: getTetV
    .. automethod:: setTetV
    .. automethod:: getTetVClamped
    .. automethod:: setTetVClamped
    
    **Triangular Data Access**
    
    .. automethod:: getTriArea
    .. automethod:: getTriSpecDefined
    .. automethod:: getTriCount
    .. automethod:: setTriCount
    .. automethod:: getTriAmount
    .. automethod:: setTriAmount
    .. automethod:: getTriClamped
    .. automethod:: setTriClamped
    .. automethod:: getTriSReacK
    .. automethod:: setTriSReacK
    .. automethod:: getTriSReacActive
    .. automethod:: setTriSReacActive
    .. automethod:: getTriSReacC
    .. automethod:: getTriSReacH
    .. automethod:: getTriSReacA
    .. automethod:: getTriSDiffD
    .. automethod:: setTriSDiffD
    .. automethod:: getTriV
    .. automethod:: setTriV
    .. automethod:: getTriVClamped
    .. automethod:: setTriVClamped
    .. automethod:: getTriOhmicI
    .. automethod:: getTriGHKI
    .. automethod:: getTriI
    .. automethod:: setTriIClamp
    .. automethod:: getTriVDepSReacActive
    .. automethod:: setTriVDepSReacActive

    **Vertex Data Access**
    
    .. automethod:: getVertV
    .. automethod:: setVertV
    .. automethod:: getVertVClamped
    .. automethod:: setVertVClamped
    .. automethod:: setVertIClamp

    **Membrane Data Access**
    
    .. automethod:: setMembPotential
    .. automethod:: setMembCapac
    .. automethod:: setMembVolRes
    .. automethod:: setMembRes
    
    **Batch Data Access**
    
    .. automethod:: getBatchTetCounts
    .. automethod:: getBatchTriCounts
    .. automethod:: getBatchTetCountsNP
    .. automethod:: getBatchTriCountsNP
    .. automethod:: sumBatchTetCountsNP
    .. automethod:: sumBatchTriCountsNP
    .. automethod:: sumBatchTriGHKIsNP
    .. automethod:: sumBatchTriOhmicIsNP

    **Region of Interest functions**
    
    .. automethod:: getROITetCounts
    .. automethod:: getROITriCounts
    .. automethod:: getROITetCountsNP
    .. automethod:: getROITriCountsNP
    .. automethod:: getROIVol
    .. automethod:: getROIArea
    .. automethod:: getROICount
    .. automethod:: setROICount
    .. automethod:: getROIAmount
    .. automethod:: getROIConc
    .. automethod:: setROIClamped
    .. automethod:: setROIReacK
    .. automethod:: setROISReacK
    .. automethod:: setROIDiffD
    .. automethod:: setROIReacActive
    .. automethod:: setROISReacActive
    .. automethod:: setROIDiffActive
    .. automethod:: setROIVDepSReacActive
    .. automethod:: getROIReacExtent
    .. automethod:: resetROIReacExtent
    .. automethod:: getROISReacExtent
    .. automethod:: resetROISReacExtent
    .. automethod:: getROIDiffExtent
    .. automethod:: resetROIDiffExtent

