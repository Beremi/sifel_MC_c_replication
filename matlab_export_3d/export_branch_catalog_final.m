function export_branch_catalog_final(outdir)

  if (nargin < 1) || isempty(outdir)
    outdir = fullfile(fileparts(pwd), 'replication_output_3d');
  end

  matlab_dir = fullfile(outdir, 'matlab');
  ensure_directory(outdir);
  ensure_directory(matlab_dir);

  young = 20000;
  poisson = 0.49;
  cohesion = 50;
  phi = pi / 9;

  E_final = [
     1.0e-5, -5.0e-3, -5.0e-3, -5.0e-3,  0.0
     0.0,     0.0,     5.0e-3, -5.0e-3,  0.0
     0.0,     5.0e-3,  5.0e-3,  1.0e-2,  5.0e-3
     0.0,     0.0,     0.0,      0.0,     0.0
     0.0,     0.0,     0.0,      0.0,     0.0
     0.0,     0.0,     0.0,      0.0,     0.0
  ];
  Ep_prev = zeros(6, size(E_final, 2));

  [S, DS, return_type, Ep_final, eig_values, sigma_values] = ...
    constitutive_problem_3d(E_final, Ep_prev, young, poisson, cohesion, phi);

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
  metadata.shear = young / (2 * (1 + poisson));
  metadata.bulk = young / (3 * (1 - 2 * poisson));
  metadata.lame = metadata.bulk - 2 * metadata.shear / 3;
  metadata.c_bar = 2 * cohesion * cos(phi);
  metadata.sin_phi = sin(phi);
  metadata.ordering = 3;

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
