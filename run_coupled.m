% This is a Matlab routine to call the OpenAeroStruct coupled system module
% written in Python.

% Setup mesh and coupled system parameters
n_inboard = 4;  % number of inboard points
n_outboard = 6; % number of outboard points
[mesh, params] = coupled_setup(n_inboard, n_outboard);

M = 25;      % max iterations
tol = 1e-6;  % relative convergence tolerance
loads = zeros(size(mesh,1)/2,6);  % initial loads array
converged = false;

% Create anonymous functions for easier use
aero = @(mesh) coupled_aero(mesh, params);
struct = @(loads) coupled_struct(loads, params);

% Check the residual of the Froebinius norm to determine if the matrix 
% has converged.
chkcnv = @(pmat,mat) logical((norm(mat,'fro')-norm(pmat,'fro'))/...
    norm(mat,'fro') < tol); 

for i = 1:M
    prevL = loads; prevM = mesh;
    
    loads = aero(mesh); % Run Aero
    mesh = struct(loads); % Run Struct
    
    % check convergence
    if chkcnv(prevL,loads) && chkcnv(prevM,mesh)
        converged = true;
        break
    end
    
end

if converged 
    fprintf('Converged after %i iterations\n',i);
else
    fprintf('Failed to converge after %i iterations\n',i);
end

% Display Loads and Mesh solutions 
fprintf('Loads...\n')
disp(loads)
fprintf('Mesh...\n')
disp(mesh)
% fprintf('Resids...\n')

% Plot the data
coupled_plotdata(mesh,loads);