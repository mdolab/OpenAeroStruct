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

try
pyversion 'C:/users/samfriedman/repos/OpenAeroStruct/venv/Scripts/python.exe';  % use python executable in virtual environment
catch
end
% MUST RUN THIS BEFORE LOADING PYTHON !!!!
try  % this only works in Unix systems
py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND
catch
end
OAS_PATH = '../';  % relative path to main OpenAeroStruct directory
%OAS_PATH = 'C:\Users\samfriedman\repos\OpenAeroStruct';
P = py.sys.path;
if count(P,OAS_PATH) == 0
    insert(P,int32(0),OAS_PATH);
end

oas_mod = py.importlib.import_module('OAS_run');

% convert cell array of lifting surfaces to python list
wing = struct;
wing.name = 'wing';
wing.num_y = int32(7);
wing.num_x = int32(3);
wing.wing_type = 'CRM';
wing.CD0 = 0.015;
wing.symmetry = true;
wing.num_twist_cp = int32(1);
wing.num_thickness_cp = int32(1);
wing.exact_failure_constraint = true;
wing.span_cos_spacing = 0.5;
wing.des_vars = py.list({'thickness_cp','twist_cp','dihedral','sweep','span','chord_cp','taper',...
                         'xshear_cp','yshear_cp','zshear_cp','radius_cp'});

tail = struct;
tail.name = 'tail';
tail.num_y = int32(7);
tail.num_x = int32(3);
tail.span = 20.0;
tail.root_chord = 5.0;
tail.wing_type = 'rect';
tail.offset = [50, 0, 5];
tail.twist_cp = [-9.5];
tail.exact_failure_constraint = true;

surf_list = {wing, tail};

prob_struct = struct;
prob_struct.type = 'aerostruct';
prob_struct.optimize = false;
prob_struct.with_viscous = true;
prob_struct.cg = py.numpy.array([70., 0., 15.]);
prob_struct.desvars = py.list({'alpha'});
prob_struct.print_level = int64(2);
%prob_struct.record_db = false;

OASobj = OAS_setup(prob_struct, surf_list);  % call matlab wrapper
fprintf('setup complete \n');

% design variables for analysis
desvars ={...
	'alpha',-5,...
	'wing.twist_cp',-10,...
	'wing.thickness_cp',0.001,...
	'wing.dihedral',-1,...
	'wing.sweep',-1,...
	'wing.span',60,...
	'wing.taper',1.5,...
	'wing.chord_cp',0.9};
 %'tail.twist_cp',[2.3],'wing.thickness_cp',[5,4]};
%for i = 1:2:length(desvars)
%     fprintf('%s, %f \n',[desvars{i}, desvars{i+1}]);
%end
disp(desvars)
output = OAS_run(desvars,OASobj);  % call matlab wrapper

% verify design point satisfies constraints
% tol = 0.1;
% if (abs(output.L_equals_W)< tol) && (output.wing_failure < 0) && (output.wing_thickness_intersects < 0) && (output.tail_failure < 0) && (output.tail_thickness_intersects < 0)
%   meets_constraints = true;
% else
%   meets_constraints = false;
% end

disp(output)
