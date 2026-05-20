function export_loading_process_final(outdir, require_all_branches)

  if (nargin < 1) || isempty(outdir)
    outdir = fullfile(fileparts(pwd), 'replication_output_3d', 'full');
  end
  if (nargin < 2) || isempty(require_all_branches)
    require_all_branches = true;
  end

  matlab_dir = fullfile(outdir, 'matlab');
  mesh_dir = fullfile(outdir, 'mesh');
  ensure_directory(outdir);
  ensure_directory(matlab_dir);
  ensure_directory(mesh_dir);

  global young poisson cohesion phi
  global shear bulk lame c_bar sin_phi specific_weight
  global x1 x2 x3 y1 y2 z N_h
  global tolerance it_max settle_max
  global step_max d_zeta_min d_zeta_init d_settle
  global target_displacement post_plastic_d_zeta
  global nv nt nint nstrain coord elem surf Dir WEIGHT

  result = loading_process;
  last_success = result.last_success;
  response = last_success.response;

  if ~any(response.return_type > 0)
    error('The final 3D full-path export has no plastic integration points.');
  end

  metadata = struct();
  metadata.nt = nint;
  metadata.nv = nv;
  metadata.nelem = nt;
  metadata.nstrain = nstrain;
  metadata.accepted_steps = last_success.accepted_steps;
  metadata.final_printed_step = last_success.printed_step;
  metadata.final_zeta = last_success.zeta;
  metadata.final_settlement = last_success.settlement;
  metadata.final_max_displacement = last_success.max_displacement;
  metadata.last_iter = last_success.iter;
  metadata.final_d_zeta = last_success.d_zeta_after_accept;
  metadata.exit_reason_code = result.exit_reason_code;
  metadata.young = young;
  metadata.poisson = poisson;
  metadata.cohesion = cohesion;
  metadata.phi = phi;
  metadata.shear = shear;
  metadata.bulk = bulk;
  metadata.lame = lame;
  metadata.c_bar = c_bar;
  metadata.sin_phi = sin_phi;
  metadata.specific_weight = specific_weight;
  metadata.x1 = x1;
  metadata.x2 = x2;
  metadata.x3 = x3;
  metadata.y1 = y1;
  metadata.y2 = y2;
  metadata.z = z;
  metadata.N_h = N_h;
  metadata.tolerance = tolerance;
  metadata.it_max = it_max;
  metadata.step_max = step_max;
  metadata.d_zeta_min = d_zeta_min;
  metadata.d_zeta_init = d_zeta_init;
  metadata.post_plastic_d_zeta = post_plastic_d_zeta;
  metadata.target_displacement = target_displacement;
  metadata.d_settle = d_settle;
  metadata.settle_max = settle_max;
  metadata.plastic_points = sum(response.return_type > 0);
  metadata.elastic_points = sum(response.return_type == 0);
  metadata.smooth_points = sum(response.return_type == 1);
  metadata.left_edge_points = sum(response.return_type == 2);
  metadata.right_edge_points = sum(response.return_type == 3);
  metadata.apex_points = sum(response.return_type == 4);

  criteria = mc_surface_coverage(response.sigma_values, response.return_type, sin_phi, c_bar);
  metadata.max_plastic_yield_residual = criteria.max_plastic_yield_residual;
  metadata.max_left_edge_gap = criteria.max_left_edge_gap;
  metadata.max_right_edge_gap = criteria.max_right_edge_gap;
  metadata.max_apex_principal_gap = criteria.max_apex_principal_gap;

  if require_all_branches
    branch_counts = [metadata.elastic_points, metadata.smooth_points, ...
                     metadata.left_edge_points, metadata.right_edge_points, ...
                     metadata.apex_points];
    if any(branch_counts == 0)
      error('The final 3D full-path export does not cover all MC return branches.');
    end
  end

  write_metadata(fullfile(outdir, 'meta.txt'), metadata);

  write_ascii_matrix(fullfile(mesh_dir, 'coord.txt'), coord);
  write_ascii_matrix(fullfile(mesh_dir, 'elem.txt'), elem);
  write_ascii_matrix(fullfile(mesh_dir, 'surf.txt'), surf);
  write_ascii_matrix(fullfile(mesh_dir, 'Dir.txt'), double(Dir));
  write_ascii_matrix(fullfile(mesh_dir, 'WEIGHT.txt'), WEIGHT);

  write_ascii_matrix(fullfile(matlab_dir, 'zeta_hist.txt'), result.zeta_hist);
  write_ascii_matrix(fullfile(matlab_dir, 'settle_hist.txt'), result.settle_hist);
  write_ascii_matrix(fullfile(matlab_dir, 'max_displacement_hist.txt'), result.max_displacement_hist);
  write_ascii_matrix(fullfile(matlab_dir, 'iter_hist.txt'), result.iter_hist);
  write_ascii_matrix(fullfile(matlab_dir, 'plastic_count_hist.txt'), result.plastic_count_hist);
  write_ascii_matrix(fullfile(matlab_dir, 'return_type_count_hist.txt'), result.return_type_count_hist);

  write_ascii_matrix(fullfile(matlab_dir, 'U_final.txt'), last_success.U);
  write_ascii_matrix(fullfile(matlab_dir, 'Ep_prev.txt'), last_success.Ep_prev);
  write_ascii_matrix(fullfile(matlab_dir, 'Ep_final_matlab.txt'), last_success.Ep_final);
  write_ascii_matrix(fullfile(matlab_dir, 'E_final.txt'), response.E_sifel);
  write_ascii_matrix(fullfile(matlab_dir, 'S_matlab.txt'), response.S_sifel);
  write_ascii_matrix(fullfile(matlab_dir, 'DS_matlab.txt'), response.DS_sifel);
  write_ascii_matrix(fullfile(matlab_dir, 'return_type_matlab.txt'), response.return_type);
  write_ascii_matrix(fullfile(matlab_dir, 'eig_matlab.txt'), response.eig_values);
  write_ascii_matrix(fullfile(matlab_dir, 'sigma_matlab.txt'), response.sigma_values);
  write_ascii_matrix(fullfile(matlab_dir, 'plastic_mask.txt'), double(response.return_type > 0));
  write_ascii_matrix(fullfile(matlab_dir, 'plastic_strain_norm.txt'), sqrt(sum(last_success.Ep_final.^2, 1)));

end

function criteria = mc_surface_coverage(sigma_values, return_type, sin_phi, c_bar)

  plastic = return_type > 0;
  left_edge = return_type == 2;
  right_edge = return_type == 3;
  apex = return_type == 4;

  criteria = struct();
  criteria.max_plastic_yield_residual = 0;
  criteria.max_left_edge_gap = 0;
  criteria.max_right_edge_gap = 0;
  criteria.max_apex_principal_gap = 0;

  if any(plastic)
    f = (1 + sin_phi) * sigma_values(1, plastic) - ...
        (1 - sin_phi) * sigma_values(3, plastic) - c_bar;
    criteria.max_plastic_yield_residual = max(abs(f));
  end

  if any(left_edge)
    criteria.max_left_edge_gap = max(abs(sigma_values(1, left_edge) - ...
                                         sigma_values(2, left_edge)));
  end

  if any(right_edge)
    criteria.max_right_edge_gap = max(abs(sigma_values(2, right_edge) - ...
                                          sigma_values(3, right_edge)));
  end

  if any(apex)
    gap_12 = abs(sigma_values(1, apex) - sigma_values(2, apex));
    gap_23 = abs(sigma_values(2, apex) - sigma_values(3, apex));
    criteria.max_apex_principal_gap = max([gap_12, gap_23]);
  end

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
  fields = fieldnames(metadata);
  for i = 1:numel(fields)
    key = fields{i};
    fprintf(fid, '%s %.17e\n', key, metadata.(key));
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
