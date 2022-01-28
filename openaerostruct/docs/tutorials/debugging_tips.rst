.. _Debugging Tips:

Debugging Tips
==============

Here are some common issues and possible solutions for them.

- When debugging, you should use the smallest mesh size possible to minimize time invested.

- Start with a simple optimization case and gradually add design variables and constraints while examining the change in optimized result. The design trade-offs in a simpler optimization problem should be easier to understand.

- Use `plot_wing` often to visualize your results. Make sure your wing geometry and structure are set up like you expect.

- After running an optimization case, always check to see if your constraints are satisfied. pyOptSparse prints out this information automatically, but if you're using Scipy's optimizer, you may need to manually print the constraints.

- Check out https://openmdao.org/newdocs/versions/latest/features/recording/case_reader_data.html for more info about how to access the saved data in the outputted `.db` file.

- If you are unsure of where a parameter or unknown is within the problem, view the N2 diagram. You can generate an html file by writing `openmdao.n2(prob)` in the runscript. You can use the search function to look for specific variables.

- Use the `thickness_intersects` constraint when varying thickness to maintain a physically possible solution.