function result = test_suite(varargin)
% Run test suite for OpenAeroStruct implementation in Matlab.
% Use optional arguments to only run certain tests. Can use multiple tags.
%   'aerostruct' = run aerostructural tests
%   'aero'       = run aero tests
%   'struct'     = run structural tests
%
% Example: test_suite('aerostruct','aero')

import matlab.unittest.TestSuite
import matlab.unittest.selectors.HasTag

tags = {};
if nargin > 0
    % only run tests with the tags included in the optional input arguments
    if any(strcmpi(varargin,'aerostruct'))
        tags{end+1} = 'Aerostruct';
    end
    if any(strcmpi(varargin,'aero'))
        tags{end+1} = 'Aero';
    end
    if any(strcmpi(varargin,'struct'))
        tags{end+1} = 'Struct';
    end
    if any(strcmpi(varargin,'misc'))
        tags{end+1} = 'Misc';
    end
    if any(strcmpi(varargin,'other'))
        tags{end+1} = 'Other';
    end
end

% Load Python
pyversion;

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

% set fortran flag
fortran_flag = py.OpenAeroStruct.run_classes.fortran_flag;

suite = TestSuite.fromClass(?Test_OpenAeroStruct);

% only select certain tests if user supplies optional argument
if ~isempty(tags)
    num_tags = length(tags);
    tests = matlab.unittest.Test.empty;
    for i = 1:num_tags
        tests = [tests, suite.selectIf(HasTag(tags{i}))]; %#ok<AGROW>
    end
else
    tests = suite;
end

% Remove Fortran tests if Fortran library not compiled
if ~fortran_flag
    tests = tests.selectIf(~HasTag('Fortran'));
end

% Run the tests, display any tests that fail or error
result = run(tests);

end