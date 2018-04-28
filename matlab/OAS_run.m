function output = OAS_run(desvars, OASobj)

%fprintf('in OAS_run \n');
%class(desvars)
%disp(desvars);
%fprintf('OASobj:\n');
%disp(OASobj);
%disp(OASobj.desvars);
%disp(OASobj.prob);
% Convert matlab cell array of design variables into python dict
py_desvars = py.dict;
for i = 1:2:length(desvars)    % desvars must have even length
    update(py_desvars, py.dict(pyargs( desvars{i}, desvars{i+1} ) ) );
end

% run OAS_run.OAS_run function and return python dict output
%class(py_desvars)
%disp(py_desvars)
py_output = py.OAS_run.OAS_run_matlab(py_desvars,OASobj);

% convert python dict to matlab struct
output = struct(py_output);
% convert values of struct array from python objects to matlab
fnames = fieldnames(output);
%save('output');
for idx = 1:length(fnames)
    fname_str = fnames{idx};
    % replace '.' in field name with '_' to work with matlab struct object
    fname_str = strrep(fname_str,'.','_');
    % fprintf('idx: %d fname{idx}: %s \n',[idx,fnames{idx}]);
    val = np2mat(output.(fname_str));
    % fprintf([fnames{idx},'=%f \n'],val)
    % replace constraint arrays with maximum values
    if any(strfind(fname_str,'failure')) || any(strfind(fname_str,'thickness_intersects'))
	output.(fname_str) = max(val(:));
    else
        output.(fname_str) = val;
    end
end

end
