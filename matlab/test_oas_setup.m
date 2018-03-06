% TEST TEST TEST

% Matlab routine to call OAS_run.OAS_setup

% This function sets up the coupled aerostructural system from the
% OpenAeroStruct analysis/optimization problem in Python. See documentation 
% and source code at https://www.github.com/samtx/OpenAeroStruct
% INPUTS:
%  prob_dict  =  (struct) Specifies design variables and problem
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
py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND

OAS_PATH = '/general/home/samfriedman/OpenAeroStruct';
%OAS_PATH = 'C:\Users\samfriedman\repos\OpenAeroStruct';
P = py.sys.path;
if count(P,OAS_PATH) == 0
    insert(P,int32(0),OAS_PATH);
end

py.importlib.import_module('OAS_run');

% convert cell array of lifting surfaces to python list
wing = struct;
wing.name = 'wing';
wing.num_y = int32(7);
wing.num_x = int32(2);
wing.wing_type = 'CRM';
wing.CD0 = 0.015;
wing.symmetry = true;
wing.num_twist_cp = int32(2);
wing.num_thickness_cp = int32(2);

tail = struct;
tail.name = 'tail';
tail.num_y = int32(7);
tail.num_x = int32(2);
tail.span = 20.0;
tail.root_chord = 5.0;
tail.wing_type = 'rect';
tail.offset = [50, 0, 5];
tail.twist_cp = [-9.5];

surf_list = {wing, tail};

py_surf_list = py.list();
for i = 1:length(surf_list)
    py_surf_list.append(surf_list{i});  % append struct (dict) to python list
end

prob_struct = struct;
prob_struct.type = 'aerostruct';
prob_struct.optimize = false;
prob_struct.with_viscous = true;
prob_struct.cg = py.numpy.array([30., 0., 5.]);
prob_struct.desvars = py.list


OASobj = py.OAS_run.OAS_setup(prob_struct, py_surf_list); 


% design variables for analysis
desvars ={'alpha',3.2,}; %'tail.twist_cp',[2.3],'wing.thickness_cp',[5,4]};


py_desvars = py.dict;
for i = 1:2:length(desvars)    % desvars must have even length
    update(py_desvars, py.dict(pyargs( desvars{i}, desvars{i+1} ) ) );
end

output = py.OAS_run.OAS_run(py_desvars,OASobj);

disp(output)















% dm = np2mat(def_mesh); % convert numpy ndarray to matlab array
