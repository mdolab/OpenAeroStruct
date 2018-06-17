# OpenAeroStruct Matlab Implementation


## Requirements
-  Matlab 2014b or newer

## Setup
Calling OpenAeroStruct from a Matlab script or function requires a few steps in succession:

1. If you are using a Unix or Linux system, then you must additionally set a library variable before loading the Python interpreter into Matlab. This variable affects how Python opens shared object libraries in Unix. This step can be made environment-agnostic with a `try/catch/end` statement.
```matlab
try
	  % On Unix/Linux systems this setting is required otherwise Matlab crashes
	  py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND
catch
end
```
2. Load the Python interpreter using `pyversion`. If you are using a virtual environment to contain the dependencies of OpenAeroStruct, then you need to specify the virtual copy of the Python executable. If you are using [virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/#lower-level-virtualenv) on Linux stored in a directory named `venv`, then the command is `pyversion <path-to-code>/OpenAeroStruct/venv/bin/python`
3. Add the OpenAeroStruct directory to the PYTHONPATH.
You can add a check if OpenAeroStruct is on your PYTHONPATH and add it if it isn't with the statements:
```matlab
OAS_PATH = full/path/to/OpenAeroStruct/directory;
if count(py.sys.path,OAS_PATH) == 0
	  insert(py.sys.path,int32(0),OAS_PATH);
end
```

Test the OpenAeroStruct setup in Matlab by running the `test_suite.m` script in the `matlab` directory.

#### References
- [System requirements for calling Python functions from Matlab](https://www.mathworks.com/help/matlab/matlab_external/system-and-configuration-requirements.html)
- [Create Python object in Matlab](https://www.mathworks.com/help/matlab/matlab_external/create-object-from-python-class.html)
- [Set DL open flags (Unix only)](https://docs.python.org/2/library/sys.html#sys.setdlopenflags)
- [pyversion](https://www.mathworks.com/help/matlab/ref/pyversion.html?s_tid=doc_ta)
- [virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/#lower-level-virtualenv)
- [PYTHONPATH](https://docs.python.org/2/using/cmdline.html#envvar-PYTHONPATH)
<!-- Import the OpenAeroStruct module. Help for this step can be found on the Mathworks documentation [Call User Defined Custom Module](https://www.mathworks.com/help/matlab/matlab_external/call-user-defined-custom-module.html). Use the import command `py.importlib.import_module('OpenAeroStruct');`.  -->

## Usage
There are several ways to run the OpenAeroStruct models from Matlab and there are several ways to retrieve data after the model is run. See `run_vlm.m`, `run_spatialbeam.m`, and `run_aerostruct.m` for the equivalent Matlab scripts to the Python examples in the parent directory.


1. Create `prob_dict` for the OASProblem object with a Matlab `struct` variable.
```Matlab
prob_dict = struct;
prob_dict.type = 'aerostruct';
```
2. Create `OASProblem` object from `run_classes.py`.  Use the `prob_dict` struct variable from before to instantiate the OASProblem with the command
```Matlab
OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict)
```
3. Add surfaces to OAS_prob with a `struct` variable.
```Matlab
surf_dict = struct;
surf_dict.num_y = 7;
surf_dict.num_x = 2;
surf_dict.wing_type = 'CRM';
surf_dict.num_twist_cp = 2;
surf_dict.num_thickness_cp = 2;
OAS_prob.add_surface(surf_dict);
```
4. Add design variables, constraints, and objectives for an optimization problem. If you only want to run a model analysis, you only need to add design variables but don't need optional design parameters.
For optimization:
```Matlab
OAS_prob.add_desvar('alpha', pyargs('lower',-10., 'upper',10.));
OAS_prob.add_desvar('wing.twist_cp', pyargs('lower',-15.,'upper',15.));
OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.01,
OAS_prob.add_constraint('wing_perf.failure',pyargs('lower',0.01, 'upper',0.5, 'scaler',1e2));
OAS_prob.add_constraint('wing_perf.thickness_intersects', pyargs('upper',0.));
OAS_prob.add_constraint('L_equals_W', pyargs('equals', 0.));
OAS_prob.add_objective('fuelburn', pyargs('scaler', 1e-5));
'upper',0.5, 'scaler',1e2));
```
For analysis only:
```matlab
OAS_prob.add_desvar('alpha');
OAS_prob.add_desvar('wing.twist_cp');
OAS_prob.add_desvar('wing.thickness_cp');
```
5. Set initial points for design variables or inputs for analysis. Use the `OASProblem.set_var()` function to do this in Matlab.
```matlab
OAS_prob.set_var('alpha', 2.5);
OAS_prob.set_var('wing.twist_cp', [-12.0, -3.0]);
OAS_prob.set_var('wing.thickness_cp', [0.03, 0.07]);
```
6. Set up the model.
```Matlab
OAS_prob.setup();
```
6. Run the model.
```Matlab
OAS_prob.run();
```
You can also use a cell array to set the design variable values when the model is run. Use the pattern `{name1, value1, name2, value2, ...}`.
```matlab
OAS_prob.run(pyargs(...
	'alpha', 2.5,...
	'wing.twist_cp', [-12.0, -3.0],...
	'wing.thickness_cp', [0.03, 0.07]));
% alternatively
input = {'alpha', 2.5,...
	'wing.twist_cp', [-12.0, -3.0],...
	'wing.thickness_cp', [0.03, 0.07]};
OAS_prob.run(pyargs(input{:}));
```
To return a struct variable with values for design variables, constraints, and objectives. Include the keyword input `{'matlab', true}` to format the output dictionary to be converted to a Matlab struct. This converts the periods in variable names to underscores. See below for information on using `pyargs()` for keyword arguments.
```Matlab
output = struct(OAS_prob.run(pyargs('matlab', true)));
```
7. Get values for design variables, constraints, and objectives.
```Matlab
fuelburn = OAS_prob.get_var('fuelburn');
alpha = OAS_prob.get_var('alpha');
twist_cp = OAS_prob.get_var('wing.twist_cp');
thickness_cp = OAS_prob.get_var('wing.thickness_cp');
```
If using an output struct from the model run:
```Matlab
fuelburn = output.fuelburn;
alpha = output.alpha;
twist_cp = output.wing_twist_cp;
thickness_cp = output.wing_thickness_cp;
```

<!-- To run the model
1. Run optimization
2. Set design variables manually and run analysis

Python commands

To get the data
1. Use the `OASProblem` object to retrieve data
2. -->

## Passing variable data between Matlab and Python
Matlab automatically converts some native Matlab variable classes to native Python types, and some Python variables can be converted to Matlab variables with built-in Matlab functions. Matlab supports passing 1-dimensional arrays to Python and OpenAeroStruct has built-in input validation that will automatically convert the appropriate parameters from floats to integers, so no additional effort is required in Matlab for to submit variables that are  1-dimensional arrays, scalar floats, or scalar integers.

#### Multidimensional Arrays
For multidimensional arrays, the Matlab variable must be converted to a Numpy array with the utility function `mat2np()` before being passed to Python. For an example, see the how the `loads` data is added to the wing surface in `run_spatialbeam.m`.
```matlab
num_y = 11;
loads = zeros(floor((num_y+1)/2), 6);
loads(:,2) = 1e5;        % This is a Matlab double (6,6) array
loads = mat2np(loads);   % Now it is a Numpy (6,6) ndarray
```
When returning data from OpenAeroStruct, Matlab will automatically convert Python scalar floats and integers to native Matlab classes. All array variables remain as a [Numpy ndarray](https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html). The utility function `np2mat()` converts the Numpy array to a native Matlab array of double precision floats. This function works for scalar values as well as 1-dimensional and multi-dimensional arrays.
For example, to return the twist control points for the wing lifting surface for use in Matlab, enter
```Matlab
wing_twist_cp = OAS_prob.get_var('wing.twist_cp');  % This is a Numpy ndarray
wing_twist_cp = np2mat(wing_twist_cp);      % Now it is a Matlab double array
```

#### Python keyword arguments
The OpenAeroStruct functions `add_desvar()`, `add_constraint()`, `add_objective()`, and `run()` can accept optional Python keyword arguments. Matlab scripts can pass these keyword arguments to Python functions with the Matlab `pyargs()` function, which groups these keyword arguments together. The function accepts pairs of inputs in the style `(key1, value1, key2, value2, ...)`. The order of the keyword pairs does not matter. 

For example, to add `wing.thickness_cp` as a design variable with bounds `[0.001, 0.5]` and scaling `1e2`, the keyword arguments to `OASProblem.add_desvar()` are

		lower = 0.001
		upper = 0.5
		scaler = 1e2

and entered in Python with the command
```Python
OASProblem.add_desvar('wing.thickness_cp', lower=0.001, upper=0.5, scaler=1e2)
```
In Matlab, use the `pyargs()` function to group the keyword arguments together.
```Matlab
OASProblem.add_desvar('wing.thickness_cp', pyargs('lower',0.001,'upper',0.5,'scaler',1e2));
```

#### Matlab struct and Python dictionary
The `run()` method to `OASProblem` returns a Python `dict` with the values for the problem design variables, constraints, objectives, and other parameters. To use this data in Matlab, the Python dictionary needs to be converted to a Matlab struct variable. In order for the conversion to work, `run()` must be called with the keyword argument `matlab = true` so that the dictionary keys are compatibile Matlab struct fieldnames. The `dict` variable can then be converted to a `struct` with the Matlab `struct()` function. The values of the struct fields are converted to Matlab variables if they are scalar but arrays remain as Numpy `ndarrays`. These `ndarrays` have to be converted to Matlab arrays by the `np2mat()` function in order to be manipulated in locally in Matlab.
```matlab
out = OAS_prob.run(pyargs('matlab',true));  % This is a Python dict
out = struct(out);		% Now this is a Matlab struct
wing_twist_cp = out.wing_twist_cp;  % This is a Numpy ndarray
wing_twist_cp = np2mat(wing_twist_cp);  % Now this is a Matlab double array
```

#### References:
- [Numpy ndarray](https://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html)
- [Pass data to Python from Matlab](https://www.mathworks.com/help/matlab/matlab_external/handle-data-returned-from-matlab-to-python.html)
- [Pass data to Matlab from Python ](https://www.mathworks.com/help/matlab/matlab_external/pass-data-to-matlab-from-python.html)
- [pyargs: Python keyword arguments](https://www.mathworks.com/help/matlab/ref/pyargs.html?s_tid=doc_ta)
- [Convert Python dict to Matlab struct](https://www.mathworks.com/help/matlab/matlab_external/convert-python-dict-type-to-matlab-structure.html)

## Troubleshooting
Many problems can be resolved by viewing the Mathworks documentation for calling Python libraries from Matlab here: [Getting Started with Python](https://www.mathworks.com/help/matlab/getting-started-with-python.html).
