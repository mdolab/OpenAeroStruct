function npary = mat2np(mat)
% MAT2NP Convert a Matlab double array to a Python Numpy ndarray.
%    NPARY = MAT2NP(MAT) returns the equivalent Numpy ndarray NPARY from
%    the Matlab double array MAT.
% 
%    MAT can be a Matlab scalar, 1-dimensional array, or n-dimensional
%    array.
%
%    If MAT is a scalar or 1-dimensional vector, then NPARY is a flattened
%    1-dimensional ndarray. If MAT is an array with 2 or greater
%    dimensions, then NPARY is an equivalently shaped ndarray.
% 
%    Complex variables are not fully supported.
%
%    Based on help from:
%    https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double
%
%    See also NP2MAT
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
else
    % N-D array, N >= 3
    % transpose first two dimensions
    mat_t = permute(mat, length(sh):-1:1);
    % pass it to Python, then reshape to python order of array shape
    npary = py.numpy.reshape(mat_t(:)', int32(fliplr(size(mat_t))));
end
end

