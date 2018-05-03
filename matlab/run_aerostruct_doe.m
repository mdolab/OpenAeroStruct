% Similar to run_aerostruct_doe.py example script

try
    % on Unix systems this setting is required otherwise Matlab crashes
    py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND
catch
end

% load Python from virtual environment with OpenMDAO 1.7.3 installed
fprintf('Load Python... \n')
[~,~,isloaded] = pyversion;
%if ~isloaded
%   pyversion 'C:\Users\Sam\repos\OpenAeroStruct\venv\Scripts\python.exe'
%end

% import OpenAeroStruct python module
OAS_PATH = '/general/home/samfriedman';
%OAS_PATH = 'C:\Users\samfriedman\repos\OpenAeroStruct';
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
OAS_prob.add_desvar('wing.span');
OAS_prob.add_desvar('wing.chord_cp');

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
%fprintf('Add design variables and constraints... \n');
%OAS_prob.add_desvar('wing.twist_cp', pyargs('lower',-15.,'upper',15.));
%OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.01, 'upper',0.5, 'scaler',1e2));
%OAS_prob.add_constraint('wing_perf.failure', pyargs('upper',0.));
%OAS_prob.add_constraint('wing_perf.thickness_intersects', pyargs('upper',0.));
%OAS_prob.add_desvar('tail.twist_cp', pyargs('lower',-15., 'upper',15.));
%OAS_prob.add_desvar('tail.thickness_cp', pyargs('lower',0.01,'upper',0.5,'scaler',1e2));
%OAS_prob.add_constraint('tail_perf.failure', pyargs('upper',0.));
%OAS_prob.add_constraint('tail_perf.thickness_intersects', pyargs('upper',0.));

% Setup problem
fprintf('Set up the problem... \n');
OAS_prob.setup()

% Generate samples
n = 5;
wing_twist_cp = linspace(-5,5,n);  	% bounds (-5, 5)  2x desvars
wing_thickness_cp = linspace(-2,2,n);	% bounds (-5, 5)  2x desvars
wing_chord_cp = linspace(0.95,1.05,n);	% bounds (0.8, 1.2)
wing_span = linspace(50,70,n);		% bounds (50, 70)
desvars = {'wing.twist_cp','wing.thickness_cp','wing.chord_cp','wing.span'};
% Get gridded vectors
[X1, X2, X3, X4, X5, X6] = ndgrid(wing_twist_cp, wing_twist_cp, wing_thickness_cp, wing_thickness_cp, wing_chord_cp, wing_span);
nsamp = numel(X1);
X = [reshape(X1, [nsamp, 1]), reshape(X2, [nsamp, 1]), reshape(X3, [nsamp, 1]), reshape(X4, [nsamp, 1]), reshape(X5, [nsamp, 1]), reshape(X6, [nsamp, 1])];

% Initialize output vectors for QOI and constraints
fuelburn = zeros(nsamp, 1);  	% QOI
failure = zeros(nsamp, 1);  	% max(wing.failure) < 0
thickness = zeros(nsamp, 1);	% max(wing.thickness_intersects) < 0
L_equals_W = zeros(nsamp,1);    % L_equals_W == 0
tol = 0.1;  % tolerance for L_equals_W == 0 constraint

% Actually run the problem
fprintf('Run the problem...  %d samples\n',nsamp);
input = cell(length(desvars)*2 + 2, 1);
input(9:10) = {'matlab',true};
tic;
for i = 1:nsamp

    if mod(i,65*50)==0
	fprintf('\n');
    end

    if mod(i,500)==0
        fprintf('%d',i);
    elseif mod(i,50)==0
	fprintf('.');
    end

    input(1:2) = {desvars{1}, X(i,1:2)};
    input(3:4) = {desvars{2}, X(i,3:4)};
    input(5:6) = {desvars{3}, X(i,5)};
    input(7:8) = {desvars{4}, X(i,6)};

    try
        output = struct(OAS_prob.run(pyargs(input{:})));
    catch
        output.fuelburn = NaN;
        output.wing_failure = NaN;
        output.wing_thickness_intersects = NaN;
    	output.L_equals_W = NaN;
    end

    fuelburn(i) = output.fuelburn;
    failure(i) = max(output.wing_failure);
    thickness(i) = max(output.wing_thickness_intersects);
    L_equals_W(i) = abs(output.L_equals_W) - tol;
end
t = toc;

%fprintf('\nFuelburn: %.4f \n', fuelburn);
fprintf('\nTime elapsed: %.4f secs\n', t);
