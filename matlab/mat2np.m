function npary = mat2np(mat)
% convert matlab matrix to python (Numpy) ndarray 
%
% Based off of https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double
%
sh = size(mat);
if any(sh(:) == 1) || any(sh(:) == 0) 
    % 1-D vector
    npary = py.numpy.array(mat(:)').flatten();
elseif length(sh) == 2
    % 2-D array
    % transpose array
    mat_t = mat';  
    % Pass array to Python as vector, then reshape to correct size
    npary = py.numpy.reshape(mat_t(:)', int32(sh));
%     mat2 = reshape(mat,1,numel(mat));  % [1, n] vector
%     npary = py.numpy.array(mat2);
%     npary = npary.reshape(int32(fliplr(sh))).transpose();  % python ndarray
else
    % N-D array, N >= 3
    % transpose first two dimensions
    mat_t = permute(mat, length(sh):-1:1);
    % pass it to Python, then reshape to python order of array shape
    npary = py.numpy.reshape(mat_t(:)', int32(fliplr(size(mat_t))));
end
end

