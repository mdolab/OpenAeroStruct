% This is a Matlab routine to test the Python aero-struct coupled system

% Setup mesh and coupled system parameters
[mesh, params] = coupled_setup(4,6);

M = 50;      % max iterations
tol = 1e-6;  % relative convergence tolerance
loads = zeros(size(mesh,1)/2,6);
converged = false;

% Create anonymous function for easier manipulation
aero = @(mesh) coupled_aero(mesh, params);
struct = @(loads) coupled_struct(loads, params);
chkcnv = @(pmat,mat) logical((norm(mat,'fro')-norm(pmat,'fro'))/...
    norm(mat,'fro') < tol);

for i = 1:M
    prevL = loads; prevM = mesh;
    % Run Aero
    loads = aero(mesh);
    % Run Struct
    mesh = struct(loads);
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
fprintf('Loads...\n')
disp(loads)
fprintf('Mesh...\n')
disp(mesh)
% fprintf('Resids...\n')

coupled_plotdata(mesh,loads);


