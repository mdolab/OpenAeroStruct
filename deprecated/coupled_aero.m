function ld = coupled_aero(dm, params)
% This function calls the aerodynamics function from the coupled.py Python
% module in the OpenAeroStruct coupled system. See documentation and source
% code at https://www.github.com/samtx/OpenAeroStruct 
% INPUTS:
%  * dm     =  (required, Matlab array) number of spanwise inboard points for
%               one side of wing mesh.
%  * params =  (required, Python dict) number of spanwise outboard points for
%               one side of wing mesh.
%
% OUTPUTS:
%    ld     =  matlab matrix of loadspoints defining the wing mesh.
%

% Call aero() function from coupled.py Python module 
def_mesh = mat2np(dm);  % convert matlab to python
loads = py.coupled.aero(def_mesh,params);
ld = np2mat(loads);     % convert python to matlab

end