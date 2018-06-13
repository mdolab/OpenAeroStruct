% Similar to run_vlm.py example script

% Script options:
%   optimize_wing = true  --> perform optimization
%   optimize_wing = false --> run a single analysis
%   use_multiple  = true  --> use multiple lifting surfaces
%   use_multiple  = false --> use a single lifting surface
optimize_wing = false;
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
prob_dict.type = 'aero';
prob_dict.optimize = optimize_wing;
prob_dict.record_db = false;  % using sqlitedict locks a process

% Instantiate OASProblem object with problem dictionary
fprintf('Create OASProblem object with prob_dwingict... \n');
OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);

%% Create a dictionary with Matlab struct to store options about the surface
surf_dict = struct;
surf_dict.num_y = 7;
surf_dict.num_x = 2;
surf_dict.wing_type = 'rect';
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
    surf_dict.span = 3.;
    surf_dict.num_y = 7;
    surf_dict.span_cos_spacing = 0.5;
    surf_dict.offset = [5., 0., .1];
    surf_dict.symmetry = true;
    fprintf('Add tail surface to problem... \n');
    OAS_prob.add_surface(surf_dict)
end

%% Add design variables, constraints, and objective for the problem
fprintf('Add design variables and constraints... \n');
OAS_prob.add_desvar('wing.twist_cp', pyargs('lower',-10., 'upper',15.));
OAS_prob.add_desvar('wing.sweep', pyargs('lower',10., 'upper',30.));
OAS_prob.add_desvar('wing.dihedral', pyargs('lower',-10., 'upper',20.));
OAS_prob.add_desvar('wing.taper', pyargs('lower',.5, 'upper',2.));
OAS_prob.add_constraint('wing_perf.CL', pyargs('equals',0.5));
OAS_prob.add_objective('wing_perf.CD', pyargs('scaler', 1e4));
if use_multiple
    OAS_prob.add_desvar('tail.twist_cp', pyargs('lower',-10., 'upper',15.));
    OAS_prob.add_desvar('tail.sweep', pyargs('lower',10., 'upper',30.));
    OAS_prob.add_desvar('tail.dihedral', pyargs('lower',-10., 'upper',20.));
    OAS_prob.add_desvar('tail.taper', pyargs('lower',.5, 'upper',2.));
    OAS_prob.add_constraint('tail_perf.CL', pyargs('equals',0.5));
end

%% Setup problem
fprintf('Set up the problem... \n');
OAS_prob.setup();

%% Actually run the problem
fprintf('Run the problem... \n');
tic;
OAS_prob.run();
t = toc;

fprintf('\nWing CL = %.13f \n', OAS_prob.getvar('wing_perf.CL'));
fprintf(  'Wing CD = %.13f \n', OAS_prob.getvar('wing_perf.CD'));
fprintf('Time elapsed: %.6f secs\n', t);
