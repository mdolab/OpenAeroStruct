function [output, meetsConstraints] = OAS_vectorize(X, desvars_list, OASobj)

% create input cell array without variable values
input = cell(1,2*length(desvars_list));
for i = 1:length(desvars_list)
    varstr = desvars_list{i};
    input{2*i-1} = varstr; 
end

% loop over rows in X and add in values and run function
output = zeros(size(X,1),1);
meetsConstraints = zeros(size(X,1),1);
LeqW_tol = 1e-2;
for i = 1:size(X,1)
    %fprintf('Sample %d \n',i);
    for j = 1:length(desvars_list)
        input{2*j} = X(i,j);
    end
    out = OAS_run(input, OASobj);
    output(i) = out.fuelburn;  % objective variable
    % check constraints
    if (out.failure<0) && (~out.thickness_intersects<0) && (abs(out.L_equals_W) < LeqW_tol)
        meetsConstraints(i) = 1;

end

end
