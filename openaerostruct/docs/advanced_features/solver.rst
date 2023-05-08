.. _Solvers Settings:

Solver Settings
===============
Users can use any of OpenMDAO's solvers for aerostructural analysis and derivative computation.
This does not means that all solvers are guaranteed to converge.
Details about the OpenMDAO's solvers can be found in the `OpenMDAO documentation <https://openmdao.org/newdocs/versions/latest/features/building_blocks/solvers/solvers.html>`_.


Default Solvers
---------------
The default nonlinear solver for aerostructural coupling is ``NonlinearBlockGS``.
For underlying aerodynamic and structural analysis, we use Scipy's LU factorization to solve linear systems (``scipy.sparse.linal.splu`` for structural FEM, and ``scipy.linalg.lu_factor/lu_solve`` for VLM).

The default linear solver for computing derivatives is ``DirectSolver``, which uses LU factorization and back substitution.
The default settings are defined in ``openaerostruct/integration/aerostruct_groups.py``.


Recommendations
---------------
For the aerostructural nonlinear solver, we recommend the default ``NonlinearBlockGS`` solver.
This solver is suitable for the circular dependency between aerodynamics and structure.
A more detailed discussion on the solver selection can be found in `this paper <http://websites.umich.edu/~mdolaboratory/pdf/Chauhan2018a.pdf>`_.

For derivative computation, the default ``DirectSolver`` works well when the mesh is coarse.
If your mesh is large (say, `nx=5` and `ny=15` or more), you may want to consider a different solver, such as ``LinearBlockGS`` or ``PETScKrylov``.
Note that these iterative linear solvers are not guaranteed to find the solution unlike the direct solver.

A good discussion on the linear solver selection can be found `here <https://openmdao.org/newdocs/versions/latest/theory_manual/setup_linear_solvers.html>`_.


Changing Solvers in Runscript
-----------------------------
You can update both the linear and nonlinear solver settings in runscript.
This should be done after calling ``prob.setup()`` and before running ``prob.run_model()`` or ``prob.run_driver()``.

Here is an example:

.. code-block::

   # First, define the problem as usual ...
   prob = Om.Problem()
   ...
   prob.model.add_subsystem('AS_point_0', AerostructPoint(), ...)
   ...

   # setup the OpenMDAO problem. This model will be initialized with the default solvers.
   prob.setup()

   # Let's try using Newton as a nonlinear solver.
   prob.model.AS_point_0.coupled.nonlinear_solver = om.NewtonSolver(solve_subsystems=True, iprint=2, maxiter=10)
   # Set a linear solver for Newton step computation (if not specified here, the linear solver attached to `prob.model.AS_point_0.coupled` will be used).
   prob.model.AS_point_0.coupled.nonlinear_solver.linear_solver = om.LinearBlockGS(iprint=-1, maxiter=30, use_aitken=True)

   # Set a linear solver to solve (a part of) the unified derivative equations.
   prob.model.AS_point_0.coupled.linear_solver = om.LinearBlockGS(iprint=1, maxiter=30, use_aitken=True)

   # Alternatively, you could use PETSc's Krylov solver with a preconditioner.
   # Note that this solver requires PETSc/petsc4py, which does not come with OpenAeroStruct by default.
   ### prob.model.linear_solver = om.PETScKrylov(assemble_jac=True, iprint=1)
   ### prob.model.options["assembled_jac_type"] = "csc"
   ### prob.model.linear_solver.precon = om.LinearRunOnce(iprint=-1)

   # run analysis or optimization.
   prob.run_model()
