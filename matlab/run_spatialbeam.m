% Similar to run_spatialbeam.py example script

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
prob_dict.type = 'struct';
prob_dict.optimize = optimize_wing;
prob_dict.record_db = false;  % using sqlitedict locks a process

% Instantiate OASProblem object with problem dictionary
fprintf('Create OASProblem object with prob_dwingict... \n');
OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);

%% Create a dictionary with Matlab struct to store options about the surface
num_y = 11;
loads = zeros(floor((num_y+1)/2), 6);
loads(:,2) = 1e5;

surf_dict = struct;
surf_dict.num_y = num_y;
surf_dict.symmetry = true;
surf_dict.loads = mat2np(loads);
fprintf('Add wing surface to problem... \n');
OAS_prob.add_surface(surf_dict);

if use_multiple
    % Multiple lifting surfaces
    surf_dict = struct;
    surf_dict.name = 'tail';
    surf_dict.span = 3.;
    surf_dict.offset = [10., 0., 0.];
    fprintf('Add tail surface to problem... \n');
    OAS_prob.add_surface(surf_dict)
end

%% Add design variables, constraints, and objective for the problem
fprintf('Add design variables and constraints... \n');
OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.001, 'upper',0.25, 'scaler',1e2));
OAS_prob.add_constraint('wing.thickness_intersects', pyargs('upper',0.));
OAS_prob.add_constraint('wing.failure', pyargs('upper',0.));
OAS_prob.add_objective('wing.structural_weight', pyargs('scaler', 1e-3));
if use_multiple
    OAS_prob.add_desvar('tail.thickness_cp', pyargs('lower',0.001,'upper',0.25,'scaler',1e2));
    OAS_prob.add_constraint('tail.thickness_intersects', pyargs('upper',0.));   
    OAS_prob.add_constraint('tail.failure', pyargs('upper',0.));
end

%% Setup problem
fprintf('Set up the problem... \n');
OAS_prob.setup();

%% Actually run the problem
fprintf('Run the problem... \n');
tic;
OAS_prob.run();
t = toc;

weight = OAS_prob.getvar('wing.structural_weight');
fprintf('\nWing structural weight: %.9f \n', weight);
fprintf('Time elapsed: %.6f secs\n', t);
