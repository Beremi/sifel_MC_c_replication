function export_branch_catalog_final(outdir)

  if (nargin < 1) || isempty(outdir)
    outdir = fullfile(fileparts(pwd), 'replication_output_3d');
  end

  matlab_dir = fullfile(outdir, 'matlab');
  ensure_directory(outdir);
  ensure_directory(matlab_dir);

  global young poisson cohesion phi
  global shear bulk lame c_bar sin_phi
  input_data;

  E_final = [
     1.0e-5, -5.0e-3, -5.0e-3, -5.0e-3,  0.0
     0.0,     0.0,     5.0e-3, -5.0e-3,  0.0
     0.0,     5.0e-3,  5.0e-3,  1.0e-2,  5.0e-3
     0.0,     0.0,     0.0,      0.0,     0.0
     0.0,     0.0,     0.0,      0.0,     0.0
     0.0,     0.0,     0.0,      0.0,     0.0
  ];
  Ep_prev = zeros(6, size(E_final, 2));

  E_upstream = sifel_to_upstream_rows(E_final);
  Ep_prev_upstream = sifel_to_upstream_rows(Ep_prev);
  [S_upstream, DS_upstream, return_type, Ep_final_upstream, eig_values, sigma_values] = ...
    constitutive_problem(E_upstream, Ep_prev_upstream);
  S = upstream_to_sifel_rows(S_upstream);
  DS = upstream_to_sifel_tangent(DS_upstream);
  Ep_final = upstream_to_sifel_rows(Ep_final_upstream);

  expected_return_type = [0, 1, 2, 3, 4];
  if ~isequal(return_type, expected_return_type)
    error('Unexpected 3D return types. Expected [%s], got [%s].', ...
          num2str(expected_return_type), num2str(return_type));
  end

  metadata = struct();
  metadata.nt = size(E_final, 2);
  metadata.young = young;
  metadata.poisson = poisson;
  metadata.cohesion = cohesion;
  metadata.phi = phi;
  metadata.shear = shear;
  metadata.bulk = bulk;
  metadata.lame = lame;
  metadata.c_bar = c_bar;
  metadata.sin_phi = sin_phi;
  metadata.sifel_ordering_code = 3;
  metadata.upstream_ordering_code = 4;

  write_metadata(fullfile(outdir, 'meta.txt'), metadata);
  write_ascii_matrix(fullfile(matlab_dir, 'E_final.txt'), E_final);
  write_ascii_matrix(fullfile(matlab_dir, 'Ep_prev.txt'), Ep_prev);
  write_ascii_matrix(fullfile(matlab_dir, 'Ep_final_matlab.txt'), Ep_final);
  write_ascii_matrix(fullfile(matlab_dir, 'S_matlab.txt'), S);
  write_ascii_matrix(fullfile(matlab_dir, 'DS_matlab.txt'), DS);
  write_ascii_matrix(fullfile(matlab_dir, 'return_type_matlab.txt'), return_type);
  write_ascii_matrix(fullfile(matlab_dir, 'eig_matlab.txt'), eig_values);
  write_ascii_matrix(fullfile(matlab_dir, 'sigma_matlab.txt'), sigma_values);

  save(fullfile(outdir, 'snapshot.mat'), ...
       'metadata', 'E_final', 'Ep_prev', 'Ep_final', 'S', 'DS', ...
       'return_type', 'eig_values', 'sigma_values');

end

function A_upstream = sifel_to_upstream_rows(A_sifel)
  sifel_to_upstream = [1, 2, 3, 6, 4, 5];
  A_upstream = A_sifel(sifel_to_upstream, :);
end

function A_sifel = upstream_to_sifel_rows(A_upstream)
  upstream_to_sifel = [1, 2, 3, 5, 6, 4];
  A_sifel = A_upstream(upstream_to_sifel, :);
end

function DS_sifel = upstream_to_sifel_tangent(DS_upstream)
  upstream_to_sifel = [1, 2, 3, 5, 6, 4];
  D_upstream = reshape(DS_upstream, 6, 6, []);
  D_sifel = D_upstream(upstream_to_sifel, upstream_to_sifel, :);
  DS_sifel = reshape(D_sifel, 36, []);
end

function ensure_directory(path)
  if exist(path, 'dir') ~= 7
    mkdir(path);
  end
end

function write_metadata(path, metadata)

  fid = fopen(path, 'w');
  if fid < 0
    error('Cannot open metadata file: %s', path);
  end

  cleanup = onCleanup(@() fclose(fid));
  names = fieldnames(metadata);
  for i = 1:numel(names)
    value = metadata.(names{i});
    if isnumeric(value)
      fprintf(fid, '%s %.17g\n', names{i}, value);
    end
  end

end

function write_ascii_matrix(path, data)

  fid = fopen(path, 'w');
  if fid < 0
    error('Cannot open output file: %s', path);
  end

  cleanup = onCleanup(@() fclose(fid));
  [nrows, ncols] = size(data);
  for i = 1:nrows
    for j = 1:ncols
      if j > 1
        fprintf(fid, ' ');
      end
      fprintf(fid, '%.17e', data(i, j));
    end
    fprintf(fid, '\n');
  end

end
