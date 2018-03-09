function npary = mat2np(mat)

% convert matlab matrix to python (Numpy) ndarray 
sh = size(mat);
if any(sh==[0,1])
    npary = py.numpy.array(mat).flatten();
else 
mat2 = reshape(mat,1,numel(mat));  % [1, n] vector
npary = py.numpy.array(mat2);
npary = npary.reshape(int32(fliplr(sh))).transpose();  % python ndarray
end

end