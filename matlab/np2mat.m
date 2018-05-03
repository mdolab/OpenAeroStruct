function mat = np2mat(npary)

% with help from: https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double

% convert python (Numpy) ndarray to matlab
try
  % if scalar
  mat = double(npary);
  % disp('A');
catch
  % if array
  sh = cellfun(@int64,cell(npary.shape));
  if length(sh) == 1
      % if a numpy 1D flattened array
      mat = double(py.array.array('d',py.numpy.nditer(npary)));
  else if length(sh) == 2
      % if a numpy 2D array
      npary2 = double(py.array.array('d',py.numpy.nditer(npary)));
      mat = reshape(npary2,fliplr(sh))';  % matlab 2d array
  else if length(sh) > 2
     % if a numpy 3D or higher dimension array
     npary3 = double(py.array.array('d',py.numpy.nditer(npary, pyargs('order','C'))));
     % fprintf('sh=%f\n',sh);
     mat = reshape(npary3, fliplr(sh));
     mat = permute(mat, [length(sh):-1:1]);  % matlab Nd array, N >= 3
  end
end
end
end
