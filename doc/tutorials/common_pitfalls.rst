.. _Common Pitfalls:

Common Pitfalls
===============

Here are some common issues and possible solutions for them.

- When debugging, you should use the smallest mesh size possible to minimize time invested.

- Start with a simple optimization case and gradually add design variables and constraints while examining the change in optimized result. The design trade-offs in a simpler optimization problem should be easier to understand.

- If you are running optimizations, it may be less computationally expensive to use `force_fd = True` within the problem dictionary.

- Use `plot_all.py` and `OptView.py` often to visualize your results. Make sure your model is set up like you expect.

- After running an optimization case, always check to see if your constraints are satisfied. pyOptSparse prints out this information automatically, but if you're using Scipy's optimizer, you may need to manually print the constraints.

- Check out http://openmdao.readthedocs.io/en/1.7.3/usr-guide/tutorials/recording.html#the-casereader for more info about how to access the saved data in the outputted `.db` file.

MDO Course Homework Tips
========================

Miscellaneous
-------------

- If you are unsure of where a parameter or unknown is within the problem, view the `.html` files that are produced when running a case. You can use the search function to look for specific variables.

- Some of the parameters may be promoted above the component level into the group level. For example, with an aerostructural case you would access the loads through `prob['coupled.wing_loads.loads']`.

- Download and use the most recent version of the OpenAeroStruct code as there may have been bug fixes for certain issues.

- Make sure all of your Python packages are up-to-date, including (but not limited to) Numpy, Scipy, Matplotlib,

Problem 6.1
-----------
- Define your aircraft using the surface and problem dictionaries, as shown in the example run scripts. In your run script, make sure to use values that are reasonable for your aircraft at the correct flight conditions. The defaults for these can be found in the `run_classes.py` file in the `get_default_prob_dict` and `get_default_surf_dict` methods (these are for a B777-sized wing). All units are SI.

- If you change an option or parameter in either the surface or problem dictionary and the results do not seem to change, check to make sure that you put the option in the correct dictionary. For example, you set the optimizer in the problem dictionary and set the CD0 in the surface dictionary.

- To define your own wing shape using OpenAeroStruct, use `wing_type : 'rect'` and then modify the wing using the geometric design variables.

- The `S_ref` variable in `surf_dict` only needs to be defined in some special cases, like aero-only planform optimization. Generally do not specify an `S_ref` value.

- You discretize your mesh using however many nodes you want. In general, `num_x = 2` and `num_y = 7` is a good starting point. You rarely need to increase `num_x` because the number of chordwise panels should not greatly influence the analysis results.

- Consult `Aircraft Design: A Conceptual Approach` by Daniel Raymer for more details about how to estimate some aircraft parameters.

- For the weights of the fuselage and other components: section 15.2.

- For CD0 estimate of fuselage and tail surfaces: sections 12.5.3 and 12.5.4.

Problem 6.2
-----------
- Make sure you choose a realistic cruise CL value based on your aircraft.

- The computed drag for the wing includes the lift-induced drag and a viscous drag estimation from empirical formulas. Your `CD0` value should include all other drag terms.

Problem 6.3
-----------
- The `structural_weight` parameter within the structural discipline is in Newtons, not kg.

- Use the `thickness_intersects` constraint when varying thickness to maintain a physically possible solution.

Problem 6.4
-----------
- To change the linear and nonlinear solvers for the aerostructural problem, look for the `ln_solver` and `nl_solver` methods within `setup_aerostruct` in `run_classes.py`. See http://openmdao.readthedocs.io/en/1.7.3/srcdocs/packages/openmdao.solvers.html for info on different solvers. Try using Newton and NLGaussSeidel on the coupled system. Check out MDO_Intro_OpenMDAO_OpenAeroStruct.pdf on Canvas for detailed information.

- Each solver has different default tolerance values. To make a valid comparison between solvers, consistently set the tolerance criteria for each solver. See the OpenMDAO documentation for the specific keywords and defaults.

Problem 6.5
-----------
- MDF is the default mode of optimization in OpenAeroStruct, so you do not need to modify anything to solve the problem using MDF.

- To create your own objective function, examine the `FunctionalBreguetRange` component within `functionals.py`. This component computes the fuel burn of the aircraft based on the calculated CL, CD, and structural weight values, along with other provided parameters. You can modify this component to also output a weighted objective between the fuel burn and the weight.

- The `FunctionalBreguetRange` component within `functionals.py` also already has the `weighted_obj` calculation inside the function. To use this objective instead, change the objective in your run script to `weighted_obj` and specify a `beta` value in the `prob_dict`. Through this method, you don't have to create your own objective function.

- See the MDO course notes for an explanation of the process for sequential optimization. The correct design variables, constraints, and objective functions are detailed there.

Problem 6.6
-----------
- To set up a structural optimization case with a certain array of loads from an aerodynamic analysis, first perform aerostructural analysis and save the outputted loads using `numpy.save()`. Then you can load the saved array using `numpy.load()` and set the loads when you initialize your structural optimization problem.

- To get the loads or displacements from an aerodynamic or structural optimization respectively, run an aerostructural analysis at that design point to perform the transfer, then input these loads or displacements into the next analysis.

- After running a simulation, the displacements can be found in `OAS_prob.prob['coupled.wing.disp']` and the loads can be found in `OAS_prob.prob['coupled.wing_loads.loads']`. Access these and save them to use in the sequential problem.

- Set these values for the loads or displacements in the `surf_dict` that you are using for your problem.
