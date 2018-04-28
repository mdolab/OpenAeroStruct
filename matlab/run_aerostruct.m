

py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND

% load Python from virtual environment with OpenMDAO 1.7.3 installed
%[~,~,isloaded] = pyversion;
%if ~isloaded
%    pyversion 'C:\Users\Sam\repos\OpenAeroStruct\oasvenv\Scripts\python.exe'
%end

% import OpenAeroStruct python module
OAS_PATH = '/general/home/samfriedman';
%OAS_PATH = 'C:\Users\samfriedman\repos\OpenAeroStruct';
P = py.sys.path;
if count(P,OAS_PATH) == 0
    insert(P,int32(0),OAS_PATH);
end
py.importlib.import_module('OpenAeroStruct');

prob_dict = struct;
prob_dict.type = 'aerostruct';
prob_dict.with_viscous = true;
prob_dict.cg = mat2np([30., 0., 5.]);
prob_dict.optimize = false;
prob_dict.record_db = false;  % using sqlitedict locks a process
prob_dict.print_level = 0;

% Instantiate problem and add default surface
OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);

% Create a dictionary to store options about the surface
surf_dict = struct;
surf_dict.name = 'wing';
surf_dict.num_y = int32(7);
surf_dict.num_x = int32(2);
surf_dict.wing_type = 'CRM';
surf_dict.CD0 = 0.015;
surf_dict.symmetry = true;
surf_dict.num_twist_cp = int32(2);
surf_dict.num_thickness_cp = int32(2);
 
% Add the specified wing surface to the problem
OAS_prob.add_surface(surf_dict);
 
%% Add design variables, constraint, and objective on the problem
%OAS_prob.add_desvar('alpha', pyargs('lower',-10., 'upper',10.));
%OAS_prob.add_constraint('L_equals_W', pyargs('equals', 0.));
%OAS_prob.add_objective('fuelburn', pyargs('scaler', 1e-5));
 
% Multiple lifting surfaces
surf_dict = struct;
surf_dict.name = 'tail';
surf_dict.num_y = int32(7);
surf_dict.num_x = int32(2);
surf_dict.span = 20.;
surf_dict.root_chord = 5.;
surf_dict.wing_type = 'rect';
surf_dict.offset = mat2np([50., 0., 5.]);
surf_dict.twist_cp = mat2np([-9.5]); %#ok<NBRAK>
% 
OAS_prob.add_surface(surf_dict)
  
% Add design variables and constraints for both the wing and tail
%OAS_prob.add_desvar('wing.twist_cp', pyargs('lower',-15.,'upper',15.));
%OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.01, 'upper',0.5, 'scaler',1e2));
%OAS_prob.add_constraint('wing_perf.failure', pyargs('upper',0.));
%OAS_prob.add_constraint('wing_perf.thickness_intersects', pyargs('upper',0.));
%OAS_prob.add_desvar('tail.twist_cp', pyargs('lower',-15., 'upper',15.));
%OAS_prob.add_desvar('tail.thickness_cp', pyargs('lower',0.01,'upper',0.5,'scaler',1e2));
%OAS_prob.add_constraint('tail_perf.failure', pyargs('upper',0.));
%OAS_prob.add_constraint('tail_perf.thickness_intersects', pyargs('upper',0.));

% Setup problem
OAS_prob.setup()

% Actually run the problem
alpha = linspace(-1,10,30);
fuelburn = zeros(size(alpha));
tic;
for i = 1:length(alpha)
    OAS_prob.setvar('alpha',alpha(i));
    
    if mod(i,10)==0
        fprintf('%i', i);
    else
        fprintf('.');
    end
    
    OAS_prob.run()
    fuelburn(i) = OAS_prob.getvar('fuelburn');
end
t = toc;

fprintf('\nFuelburn: %.4f \n', OAS_prob.getvar('fuelburn'));
fprintf('Time elapsed: %.4f secs\n', t);
plot(alpha,fuelburn);
xlabel('alpha');
ylabel('fuelburn');
