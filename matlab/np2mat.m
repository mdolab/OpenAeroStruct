function mat = np2mat(npary)
% NP2MAT Convert a Numpy ndarray to a native Matlab double array
%    MAT = NP2MAT(NPARY) returns the equivalent Matlab array MAT from the
%    Python Numpy ndarray NPARY.
%
%    NPARY can be a Python scalar, Numpy ndarray scalar, flattened Numpy
%    1-dimensional ndarray, or a Numpy n-dimensional ndarray.
%
%    Complex variables are not supported.
%
%    Based on help from:
%    https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double
%
%    See also MAT2NP
try
    % if scalar float or int value
    mat = double(npary);
catch
    % Otherwise assume it is a ndarray
    sh = cellfun(@int64,cell(npary.shape));
    if any(numel(sh) == [0,1])
        % if a scalar or numpy 1D flattened array
        mat = double(py.array.array('d',py.numpy.nditer(npary)));
    elseif length(sh) == 2
        % if a numpy 2D array
        npary2 = double(py.array.array('d',py.numpy.nditer(npary)));
        mat = reshape(npary2,fliplr(sh))';  % matlab 2d array
    elseif length(sh) > 2
        % if a numpy 3D or higher dimension array
        npary3 = double(py.array.array('d',py.numpy.nditer(npary, pyargs('order','C'))));
        mat = reshape(npary3, fliplr(sh));
        mat = permute(mat, length(sh):-1:1);  % matlab Nd array, N >= 3
    end
end
end
