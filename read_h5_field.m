function [field] = read_from_h5(file, name)
    xyz = 'xyz';
    for k = 1 : 3
        field{k} = hdf5read(file, ['/', name, '_', xyz(k), '_real']) + ...
                i * hdf5read(file, ['/', name, '_', xyz(k), '_imag']);
        field{k} = permute(field{k}, [ndims(field{k}):-1:1]);
    end
