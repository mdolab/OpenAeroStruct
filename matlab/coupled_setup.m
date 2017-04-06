function [dm, params] = coupled_setup(n_inb, n_outb, varargin)
% This function sets up the coupled aero-structural system from the
% OpenAeroStruct optimization problem in Python. See documentation and
% source code at https://www.github.com/samtx/OpenAeroStruct
% INPUTS:
%  * n_inb  =  (required, integer) number of spanwise inboard points for
%               one side of wing mesh.
%  * n_outb =  (required, integer) number of spanwise outboard points for
%               one side of wing mesh.
%    check  =  (optional, boolean) should the output of the OpenMDAO setup
%               process be printed. Default is false.
%    fname  =  (optional, string) The filename or output method that the
%               setup check should be written to. Default is sys.stdout,
%               which prints to screen.
% OUTPUTS:
%    dm     =  matlab matrix of points defining the wing mesh.
%    params =  python dict object of parameters for the coupled_aero and
%              coupled_struct functions
%

% ------------------------   TROUBLESHOOTING   -------------------------
% Python Module Not on Python Search Path
% If command is a valid Python command, make sure the Python module is on the Python search path. To test if module mymod is on the path, type:
% 
% py.importlib.import_module('mymod')
% If Python cannot find the module, MATLAB displays a Python error message.
% 
% To add mymod, in folder modpath, to the path, type:
% 
% P = py.sys.path;
% if count(P,'modpath') == 0
%     insert(P,int32(0),'modpath');
% end


nvarargs = length(varargin);
if nvarargs == 1
    check = logical(varargin{1});
    outstream = py.sys.stdout;  % print to screen
elseif nvarargs == 2
    check = logical(varargin{1});
    outstream = varargin{2};
else
    check = false;
    outstream = py.sys.stdout;
end

% Call setup() function from coupled.py Python module 
out = py.coupled.setup(n_inb, n_outb, check, outstream); 
def_mesh = out{1};   % initial mesh for wing
params = out{2};     % parameters for aero and struct modules
dm = np2mat(def_mesh); % convert numpy ndarray to matlab array

end