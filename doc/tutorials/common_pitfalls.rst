.. _Common Pitfalls:

Common Pitfalls
===============

Here are some common issues and possible solutions for them.

- When debugging, you should use the smallest mesh size possible to minimize time invested.

- Start with a simple optimization case and gradually add design variables and constraints while examining the change in optimized result. The design trade-offs in a simpler optimization problem should be easier to understand.

- If you are running optimizations, it may be less computationally expensive to use `force_fd = True` within the problem dictionary.

- Use `plot_all.py` and `OptView.py` often to visualize your results. Make sure your model is set up like you expect.

- After running an optimization case, always check to see if your constraints are satisfied. pyOptSparse prints out this information automatically, but if you're using Scipy's optimizer, you may need to manually print the constraints.

- Check out http://openmdao.readthedocs.io/en/latest/usr-guide/tutorials/recording.html#the-casereader for more info about how to access the saved data in the outputted `.db` file.

MDO Course Homework Tips
========================

.. note::
  Make sure to change the default surface and problem dictionaries to use values that are reasonable for your aircraft at the correct flight conditions. All units are SI.

- If you are unsure of where a parameter or unknown is within the problem, view the `.html` files that are produced when running a case. You can use the search function to look for specific variables.

- Download and use the most recent version of the OpenAeroStruct code as there may have been bug fixes for certain issues.

- Make sure all of your Python packages are up-to-date, including (but not limited to) Numpy, Scipy, Matplotlib,

- Use the `thickness_intersects` constraint when varying thickness to maintain a physically possible solution.

- To set up a structural optimization case with a certain array of loads from an aerodynamic analysis, first perform aerostructural analysis and save the outputted loads using `numpy.save()`. Then you can load the saved array using `numpy.load()` and set the loads when you initialize your structural optimization problem.

- To create your own objective function, examine the `FunctionalBreguetRange` component within `functionals.py`. This component computes the fuel burn of the aircraft based on the calculated CL, CD, and structural weight values, along with other provided parameters. You can modify this component to also output a weighted objective between the fuel burn and the weight.

- Change `prob_dict['g']` to alter the force of gravity for different load cases.

- The `structural_weight` parameter within the structural discipline is in Newtons, not kg.

- To change the linear and nonlinear solvers for the aerostructural problem, look for the `ln_solver` and `nl_solver` methods within `setup_aerostruct` in `run_classes.py`. See http://openmdao.readthedocs.io/en/latest/srcdocs/packages/openmdao.solvers.html for info on different solvers. Try using Newton and NLGaussSeidel on the coupled system.

.. note::
  Each solver has different default tolerance values. To make a valid comparison between solvers, consistently set the tolerance criteria for each solver. See the OpenMDAO documentation for the specific keywords and defaults.
