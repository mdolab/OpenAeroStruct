function OASobj = OAS_setup(prob_struct, surf_list)
% This function sets up the coupled aerostructural system from the
% OpenAeroStruct analysis/optimization problem in Python. See documentation 
% and source code at https://www.github.com/samtx/OpenAeroStruct
% INPUTS:
%  prob_struct  =  (struct) Specifies design variables and problem
%                parameters.
%  surf_list  =  (cell array) Cell array of structs defining the 
%                lifting surfaces. 
% OUTPUTS:
%  OASobj     =  python object of OASProblem() to pass to OAS_run.m
%
% ------------------------   TROUBLESHOOTING   -------------------------
% Python Module Not on Python Search Path
% If command is a valid Python command, make sure the Python module is on the 
% Python search path. To test if module OpenAeroStruct is on the path, type:
% 
% py.importlib.import_module('OpenAeroStruct')
% If Python cannot find the module, MATLAB displays a Python error message.
% 
% To add OpenAeroStruct, in folder OAS_PATH, to the PYTHONPATH, type:
% 
% P = py.sys.path;
% if count(P,'OAS_PATH') == 0
%     insert(P,int32(0),'OAS_PATH');
% end
%
% -----------------------------------------------------------------------

% Call OAS_setup() function from OAS_run.py Python module
% can generate correct path string by changing to directory and
% using pwd command in Matlab Command Window.


% MUST RUN THIS BEFORE LOADING PYTHON !!!!
% We might not need this if we use Intel-compiled Python/Scipy/Numpy with MKL instead of standard GCC and BLAS
% see link for more info: https://www.mathworks.com/matlabcentral/answers/327193-calling-python-module-from-matlab-causes-segmentation-fault-in-h5py
try  %only works in unix systems
py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND
catch
end
py.importlib.import_module('OAS_run');

% convert matlab cell array of structs into python list of dicts
py_surf_list = py.list();
for i = 1:length(surf_list)
    py_surf_list.append(surf_list{i});  % append struct(dict) to python list
end

% Call OAS_run.OAS_setup function and return OASProblem object
OASobj = py.OAS_run.OAS_setup(prob_struct, py_surf_list); 

end