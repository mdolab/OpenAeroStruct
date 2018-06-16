% Similar to run_aerostruct.py example script

% Script options:
%   optimize_wing = true  --> perform optimization
%   optimize_wing = false --> run a single analysis
%   use_multiple  = true  --> use multiple lifting surfaces
%   use_multiple  = false --> use a single lifting surface
optimize_wing = true;
use_multiple = true;

% Load Python
fprintf('Load Python... ')
[~,~,isloaded] = pyversion;
if isloaded
   fprintf('Python loaded.\n')
end

try
    % On Unix/Linux systems this setting is required otherwise Matlab crashes
    py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND
catch
end

% Path to OpenAeroStruct directory
OAS_PATH = py.os.path.abspath('../..');

% Add OpenAeroStruct directory to PYTHONPATH
P = py.sys.path;
if count(P,OAS_PATH) == 0
    insert(P,int64(0),OAS_PATH);
end

%% Create problem dictionary using Matlab struct object
prob_dict = struct;
prob_dict.type = 'aerostruct';
prob_dict.with_viscous = true;
prob_dict.optimize = optimize_wing;
prob_dict.record_db = false;  % using sqlitedict locks a process
prob_dict.cg = [30., 0. 5.];

% Instantiate OASProblem object with problem dictionary
fprintf('Create OASProblem object with prob_dict... \n');
OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);

%% Create a dictionary with Matlab struct to store options about the surface
surf_dict = struct;
surf_dict.num_y = 7;
surf_dict.num_x = 2;
surf_dict.wing_type = 'CRM';
surf_dict.CD0 = 0.015;
surf_dict.symmetry = true;
surf_dict.num_twist_cp = 2;
surf_dict.num_thickness_cp = 2;
fprintf('Add wing surface to problem... \n');
OAS_prob.add_surface(surf_dict);

if use_multiple
    % Multiple lifting surfaces
    surf_dict = struct;
    surf_dict.name = 'tail';
    surf_dict.num_y = 7;
    surf_dict.num_x = 2;
    surf_dict.span = 20.;
    surf_dict.root_chord = 5.;
    surf_dict.wing_type = 'rect';
    surf_dict.offset = [50., 0., 5.];
    surf_dict.twist_cp = -9.5;
    surf_dict.exact_failure_constraint = true;
    fprintf('Add tail surface to problem... \n');
    OAS_prob.add_surface(surf_dict)
end

%% Add design variables, constraints, and objective for the problem
fprintf('Add design variables and constraints... \n');
OAS_prob.add_desvar('alpha', pyargs('lower',-10., 'upper',10.));
OAS_prob.add_constraint('L_equals_W', pyargs('equals', 0.));
OAS_prob.add_objective('fuelburn', pyargs('scaler', 1e-5));
OAS_prob.add_desvar('wing.twist_cp', pyargs('lower',-15.,'upper',15.));
OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.01, 'upper',0.5, 'scaler',1e2));
OAS_prob.add_constraint('wing_perf.failure', pyargs('upper',0.));
OAS_prob.add_constraint('wing_perf.thickness_intersects', pyargs('upper',0.));
if use_multiple
    OAS_prob.add_desvar('tail.twist_cp', pyargs('lower',-15., 'upper',15.));
    OAS_prob.add_desvar('tail.thickness_cp', pyargs('lower',0.01,'upper',0.5,'scaler',1e2));
    OAS_prob.add_constraint('tail_perf.failure', pyargs('upper',0.));
    OAS_prob.add_constraint('tail_perf.thickness_intersects', pyargs('upper',0.));
end

%% Setup problem
fprintf('Set up the problem... \n');
OAS_prob.setup();

%% Actually run the problem
fprintf('Run the problem... \n');
tic;
OAS_prob.run();
t = toc;

fuelburn = OAS_prob.get_var('fuelburn');
fprintf('\nFuelburn: %f \n', fuelburn);
fprintf('Time elapsed: %.4f secs\n', t);
