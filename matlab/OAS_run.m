function output = OAS_run(desvars, OASobj)

% Convert matlab cell array of design variables into python dict
py_desvars = py.dict;
for i = 1:2:length(desvars)    % desvars must have even length
    update(py_desvars, py.dict(pyargs( desvars{i}, desvars{i+1} ) ) );
end 

% run OAS_run.OAS_run function and return python dict output
py_output = py.OAS_run.OAS_run(py_desvars,OASobj);

% convert python dict to matlab struct
output = struct(py_output);
% convert values of struct array from python objects to matlab
fnames = fieldnames(output);
%save('output');
for idx = 1:length(fnames)
    val = np2mat(output.(fnames{idx}));
    %fprintf([fnames{idx},'=%f \n'],val)
    output.(fnames{idx}) = val;
end


end
