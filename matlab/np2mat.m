function mat = np2mat(npary)

% convert python (Numpy) ndarray to matlab
try  
  % if scalar
  mat = double(npary);
catch
  % if array
  sh = cellfun(@int32,cell(npary.shape));
  if length(sh) == 1
      % if a numpy 1D flattened array
      mat = double(py.array.array('d',npary));
  else
      % if a numpy 2D array
      npary2 = double(py.array.array('d',py.numpy.nditer(npary)));
      mat = reshape(npary2,fliplr(sh))';  % matlab 2d array 
  end
end