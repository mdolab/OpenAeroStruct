function dm = coupled_struct(ld, kwargs)

% Struct
loads = mat2np(ld);
def_mesh = py.coupled.struct(loads, kwargs);  % initial mesh for wing
dm = np2mat(def_mesh);
% sh = double(py.array.array('d',def_mesh.shape));
% dm1 = double(py.array.array('d',py.numpy.nditer(def_mesh)));
% dm = reshape(dm1,fliplr(sh))';  % matlab 2d array of def_mesh

return