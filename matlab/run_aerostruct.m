% Similar to run_aerostruct.py example script

try
    % on Unix systems this setting is required otherwise Matlab crashes
    py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND
catch
end

% load Python from virtual environment with OpenMDAO 1.7.3 installed
fprintf('Load Python... \n')
[~,~,isloaded] = pyversion;
if ~isloaded
   pyversion 'C:\Users\Sam\repos\OpenAeroStruct\venv\Scripts\python.exe'
end

% import OpenAeroStruct python module
% OAS_PATH = '/general/home/samfriedman';
OAS_PATH = 'C:\Users\samfriedman\repos\OpenAeroStruct';
P = py.sys.path;
if count(P,OAS_PATH) == 0
    insert(P,int64(0),OAS_PATH);
end
fprintf('Import OpenAeroStruct module... \n');
py.importlib.import_module('OpenAeroStruct');

prob_dict = struct;
prob_dict.type = 'aerostruct';
prob_dict.with_viscous = true;
prob_dict.cg = mat2np([30., 0., 5.]);
% prob_dict.cg = mat2np([30., 0., 5.]);
prob_dict.optimize = false;
prob_dict.record_db = false;  % using sqlitedict locks a process
prob_dict.print_level = 0;

% Instantiate problem and add default surface
fprintf('Create OASProblem object with prob_dict... \n');
OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);

% Create a dictionary to store options about the surface
surf_dict = struct;
surf_dict.name = 'wing';
surf_dict.num_y = int64(7);
surf_dict.num_x = int64(2);
surf_dict.wing_type = 'CRM';
surf_dict.CD0 = 0.015;
surf_dict.symmetry = true;
surf_dict.num_twist_cp = int64(2);
surf_dict.num_thickness_cp = int64(2);

% Add the specified wing surface to the problem
fprintf('Add wing surface to problem... \n');
OAS_prob.add_surface(surf_dict);

%% Add design variables, constraint, and objective on the problem
OAS_prob.add_desvar('alpha', pyargs('lower',-10., 'upper',10.));
OAS_prob.add_constraint('L_equals_W', pyargs('equals', 0.));
OAS_prob.add_objective('fuelburn', pyargs('scaler', 1e-5));

% Multiple lifting surfaces
surf_dict = struct;
surf_dict.name = 'tail';
surf_dict.num_y = int64(7);
surf_dict.num_x = int64(2);
surf_dict.span = 20.;
surf_dict.root_chord = 5.;
surf_dict.wing_type = 'rect';
surf_dict.offset = [50., 0., 5.];
surf_dict.twist_cp = -9.5;
% surf_dict.root_chord = 5.;
% surf_dict.wing_type = 'rect';
% surf_dict.offset = mat2np([50., 0., 5.]);
% surf_dict.twist_cp = mat2np(-9.5); %#ok<NBRAK>
%
fprintf('Add tail surface to problem... \n');
OAS_prob.add_surface(surf_dict)

% Add design variables and constraints for both the wing and tail
fprintf('Add design variables and constraints... \n');
OAS_prob.add_desvar('wing.twist_cp', pyargs('lower',-15.,'upper',15.));
OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.01, 'upper',0.5, 'scaler',1e2));
OAS_prob.add_constraint('wing_perf.failure', pyargs('upper',0.));
OAS_prob.add_constraint('wing_perf.thickness_intersects', pyargs('upper',0.));
OAS_prob.add_desvar('tail.twist_cp', pyargs('lower',-15., 'upper',15.));
OAS_prob.add_desvar('tail.thickness_cp', pyargs('lower',0.01,'upper',0.5,'scaler',1e2));
OAS_prob.add_constraint('tail_perf.failure', pyargs('upper',0.));
OAS_prob.add_constraint('tail_perf.thickness_intersects', pyargs('upper',0.));

% Setup problem
fprintf('Set up the problem... \n');
OAS_prob.setup()

% Actually run the problem
fprintf('Run the problem... \n');
tic;
output = struct(OAS_prob.run(pyargs('matlab',true)));
t = toc;
fuelburn = output.fuelburn;

fprintf('\nFuelburn: %.4f \n', fuelburn);
fprintf('Time elapsed: %.4f secs\n', t);
