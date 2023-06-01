.. _Multipoint Optimization with MPI:

.. note:: This feature requires MPI, which does not come with OpenAeroStruct by default.

Parallel Multipoint Optimization using MPI
==========================================

Multipoint analysis or optimization can be parallelized to reduce the runtime.
Because each flight condition (or point) is independent, it is `embarassingly parallel`, meaning that we can easily parallelize these analyses.

Here, we will parallelize the :ref:`previous multipoint aerostructural example (Q400)<Multipoint Optimization>`.
This requires a little modification to the original serial runscript.

Runscript modifications
-----------------------
We first import MPI.
If this line does not work, make sure that you have a working MPI installation.

.. code-block::
    
    from mpi4py import MPI

You may need to turn off the numpy multithreading.
This can be done by adding the following lines before importing numpy.
The name of environment variable may be different depending on the system.

.. code-block::

    import os
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

Then, let's set up the problem in the same way as the serial runscript.

.. code-block::
    
    prob = om.Problem()
    
    # Setup problem information in indep_var_comp
    ...

    # Add AerostructGeometry to the model
    ...

Next, we need to add AS_points under a ``ParallelGroup`` instead of directly under the ``prob.model``.

.. literalinclude:: ../../tests/test_multipoint_parallel.py
    :dedent: 8
    :start-after: # [rst Setup ParallelGroup (beg)]
    :end-before: # [rst Setup ParallelGroup (end)]

After establishing variable connections and setting up the driver, we define the optimization objective and constraints.
Here, we will setup the parallel derivative computations.
In this example, we have 6 functions of interest (1 objective and 5 constraints), which would require 6 linear solves for reverse-mode derivatives in series.
Among 6 functions, 4 depend only on AS_point_0, and 2 depend only on AS_point_1.
Therefore, we can form 2 pairs and perform linear solves in parallel.
We specify ``parallel_deriv_color`` to tell OpenMDAO which function's derivatives can be solved for in parallel.

.. literalinclude:: ../../tests/test_multipoint_parallel.py
    :dedent: 8
    :start-after: # [rst Parallel deriv color setup 1 (beg)]
    :end-before: # [rst Parallel deriv color setup 1 (end)]

Furthermore, we will add another dummy (nonsense) constraint to explain how parallelization works for reverse-mode derivatives.
This dummy constraint (sum of the fuel burns from AS_point_0 and AS_point_1) depends on both AS points.
In this case, the linear solves of AS_point_0 and AS_point_1 will be parallelized.

.. literalinclude:: ../../tests/test_multipoint_parallel.py
    :dedent: 8
    :start-after: # [rst Parallel deriv color setup 2 (beg)]
    :end-before: # [rst Parallel deriv color setup 2 (end)]

Finally, let's change the linear solver from default.
This step is not necessary and not directly relevant to parallelization, but the ``LinearBlockGS`` solver works better on a fine mesh than the default ``DirectSolver``.

.. literalinclude:: ../../tests/test_multipoint_parallel.py
    :dedent: 8
    :start-after: # [rst Change linear solver (beg)]
    :end-before: # [rst Change linear solver (end)]


Complete runscript
------------------

.. embed-code::
    openaerostruct.tests.test_multipoint_parallel.Test.test_multipoint_MPI

To run this example in parallel with two processors, use the following command:

.. code-block:: bash

    $ mpirun -n 2 python <name of script>.py

Solver Outputs
--------------
The stdout of the above script would look like the following.
The solver outputs help us understand how solvers are parallelized for analysis and total derivative computations.

.. code-block:: bash

    # Nonlinear solver in parallel
    ===========================
    parallel.AS_point_0.coupled
    ===========================

    ===========================
    parallel.AS_point_1.coupled
    ===========================
    NL: NLBGS 1 ; 82168.4402 1
    NL: NLBGS 1 ; 79704.5639 1
    NL: NLBGS 2 ; 63696.5109 0.775194354
    NL: NLBGS 2 ; 68552.4805 0.860082248
    NL: NLBGS 3 ; 2269.83605 0.0276241832
    NL: NLBGS 3 ; 2641.30776 0.0331387267
    NL: NLBGS 4 ; 26.8901082 0.000327255917
    NL: NLBGS 4 ; 33.4963389 0.000420256222
    NL: NLBGS 5 ; 0.20014208 2.43575367e-06
    NL: NLBGS 5 ; 0.273747809 3.43453117e-06
    NL: NLBGS 6 ; 0.000203058798 2.47125048e-09
    NL: NLBGS 6 ; 0.00033072442 4.14937871e-09
    NL: NLBGS 7 ; 3.3285346e-06 4.05086745e-11
    NL: NLBGS 7 ; 5.16564573e-06 6.48099115e-11
    NL: NLBGS 8 ; 9.30405466e-08 1.13231487e-12
    NL: NLBGS Converged
    NL: NLBGS 8 ; 1.63279302e-07 2.04855649e-12
    NL: NLBGS 9 ; 2.01457772e-09 2.5275563e-14
    NL: NLBGS Converged

    # Linear solver for "parcon1". Derivatives of AS_point_0.fuelburn and AS_point_1.L_equals_W in parallel.
    ===========================
    parallel.AS_point_0.coupled
    ===========================

    ===========================
    parallel.AS_point_1.coupled
    ===========================
    LN: LNBGS 0 ; 180.248073 1
    LN: LNBGS 0 ; 1.17638541e-05 1
    LN: LNBGS 1 ; 0.00284457871 1.57814653e-05
    LN: LNBGS 1 ; 1.124189e-06 0.0955629836
    LN: LNBGS 2 ; 1.87700622e-08 0.00159557081
    LN: LNBGS 2 ; 4.66688449e-05 2.58914529e-07
    LN: LNBGS 3 ; 1.13549461e-11 9.65240308e-07
    LN: LNBGS Converged
    LN: LNBGS 3 ; 8.18485966e-08 4.54088609e-10
    LN: LNBGS 4 ; 9.00103905e-10 4.99369503e-12
    LN: LNBGS Converged

    # Linear solver for "parcon2". Derivatives of AS_point_0.CL and AS_point_1.wing_perf.failure in parallel.
    ===========================
    parallel.AS_point_1.coupled
    ===========================

    ===========================
    parallel.AS_point_0.coupled
    ===========================
    LN: LNBGS 0 ; 334.283603 1
    LN: LNBGS 0 ; 0.00958374526 1
    LN: LNBGS 1 ; 2.032696e-05 6.08075293e-08
    LN: LNBGS 1 ; 2.02092209e-06 0.000210869762
    LN: LNBGS 2 ; 2.3346978e-06 6.98418281e-09
    LN: LNBGS 2 ; 2.90180431e-08 3.02783956e-06
    LN: LNBGS 3 ; 4.98483883e-08 1.49120052e-10
    LN: LNBGS 3 ; 8.63240127e-11 9.0073359e-09
    LN: LNBGS Converged
    LN: LNBGS 4 ; 5.58667374e-11 1.67123774e-13
    LN: LNBGS Converged

    # Linear solver for derivatives of fuel_vol_delta.fuel_vol_delta (not parallelized)
    ===========================
    parallel.AS_point_0.coupled
    ===========================
    LN: LNBGS 0 ; 0.224468335 1
    LN: LNBGS 1 ; 3.54243924e-06 1.57814653e-05
    LN: LNBGS 2 ; 5.81181131e-08 2.58914529e-07
    LN: LNBGS 3 ; 1.01928513e-10 4.54088604e-10
    LN: LNBGS 4 ; 1.12121714e-12 4.99499023e-12
    LN: LNBGS Converged

    # Linear solver for derivatives of fuel_diff (not parallelized)
    ===========================
    parallel.AS_point_0.coupled
    ===========================
    LN: LNBGS 0 ; 0.21403928 1
    LN: LNBGS 1 ; 3.37785348e-06 1.57814653e-05
    LN: LNBGS 2 ; 5.54178795e-08 2.58914529e-07
    LN: LNBGS 3 ; 9.71927996e-11 4.54088612e-10
    LN: LNBGS Converged

    # Linear solver for derivatives of fuel_sum in parallel.
    ===========================
    parallel.AS_point_0.coupled
    ===========================

    ===========================
    parallel.AS_point_1.coupled
    ===========================
    LN: LNBGS 0 ; 360.496145 1
    LN: LNBGS 0 ; 511.274568 1
    LN: LNBGS 1 ; 0.00568915741 1.57814653e-05
    LN: LNBGS 1 ; 0.00838867553 1.64073788e-05
    LN: LNBGS 2 ; 0.00013534629 2.64723299e-07
    LN: LNBGS 2 ; 9.33376897e-05 2.58914529e-07
    LN: LNBGS 3 ; 1.00754737e-07 1.97065811e-10
    LN: LNBGS 3 ; 1.63697193e-07 4.54088609e-10
    LN: LNBGS 4 ; 2.24690253e-09 4.39470819e-12
    LN: LNBGS Converged
    LN: LNBGS 4 ; 1.80020781e-09 4.99369503e-12
    LN: LNBGS Converged


Comparing Runtime
-----------------
How much speedup can we get by parallelization?
Here, we compared the runtime for the example above (but with a finer mesh of `nx=3` and `ny=61`).
In this case, we achieved decent speedup in nonlinear analysis, but not so much in derivative computation.
The actual speedup you can get depends on your problem setups, such as number of points (flight conditions) and functions of interest.

.. list-table:: Runtime for Q400 example
   :widths: 30 35 35
   :header-rows: 1

   * - Case
     - Analysis walltime [s]
     - Derivatives walltime [s]
   * - Serial
     - 1.451
     - 5.775
   * - Parallel 
     - 0.840
     - 4.983
