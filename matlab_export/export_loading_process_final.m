function export_loading_process_final(outdir)

  if (nargin < 1) || isempty(outdir)
    outdir = fullfile(fileparts(pwd), 'replication_output');
  end

  matlab_dir = fullfile(outdir, 'matlab');
  cpp_dir = fullfile(outdir, 'cpp');
  ensure_directory(outdir);
  ensure_directory(matlab_dir);
  ensure_directory(cpp_dir);

  global shear bulk lame
  global c_bar sin_phi
  global specific_weight
  global x1 x2 x3 y1 y2 h
  global tolerance it_max settle_max
  global step_max d_zeta_min d_zeta_init d_settle
  global DEV Dev VOL iota ELAST
  global Elast Inv_ELAST IDENT Ident

  global nt nv
  global coord
  global elem
  global Dir
  global area
  global EU
  global F1
  global K_elast
  global iMALI jMALI

  input_data;
  preprocessing;

  zeta_hist = zeros(1, step_max);
  settle_hist = zeros(1, step_max);
  d_zeta = d_zeta_init;
  settle_old = 0;
  zeta_old = 0;
  U_old = zeros(2, nv);
  Ep_old = zeros(4, nt);
  step = 1;

  last_success = struct( ...
      'accepted', false, ...
      'printed_step', 1, ...
      'accepted_steps', 0, ...
      'zeta', 0, ...
      'settlement', 0, ...
      'iter', 0, ...
      'd_zeta_after_accept', d_zeta, ...
      'U', zeros(2, nv), ...
      'Ep_prev', zeros(4, nt), ...
      'Ep_final', zeros(4, nt));

  exit_reason_code = 0;

  while true

      zeta = zeta_old + d_zeta;
      F = zeta * F1;
      Ep_prev_step = Ep_old;

      [U, Ep, iter] = newton(U_old, Ep_old, F);

      if iter == it_max
          d_zeta = max(d_zeta / 2, d_zeta_min);
      else
          step = step + 1;
          zeta_old = zeta;
          zeta_hist(step) = zeta;
          settle = -U(2, nv);
          settle_hist(step) = settle;

          if (settle - settle_old) > d_settle
              warning('Too large increment of the settlement.')
              d_zeta = max(d_zeta / 2, d_zeta_min);
          end

          settle_old = settle;
          U_old = U;
          Ep_old = Ep;

          last_success.accepted = true;
          last_success.printed_step = step;
          last_success.accepted_steps = step - 1;
          last_success.zeta = zeta;
          last_success.settlement = settle;
          last_success.iter = iter;
          last_success.d_zeta_after_accept = d_zeta;
          last_success.U = U;
          last_success.Ep_prev = Ep_prev_step;
          last_success.Ep_final = Ep;

          disp([' step=', num2str(step), ', settlement=', ...
                num2str(settle), ', zeta=', num2str(zeta), ', d_zeta=', ...
                num2str(d_zeta), ', iter=', num2str(iter)])
      end

      if settle_old >= settle_max
          warning('Too large settlement at the point A.')
          exit_reason_code = 1;
          break
      end

      if d_zeta == d_zeta_min
          warning('Too small load increments.')
          exit_reason_code = 2;
          break
      end

      if step >= step_max
          warning('Maximal number of steps was achieved.')
          exit_reason_code = 3;
          break
      end

  end

  if ~last_success.accepted
      error('No successful load step was completed.');
  end

  zeta_hist = zeta_hist(1:step);
  settle_hist = settle_hist(1:step);

  E_final = zeros(3, nt);
  E_final(:) = EU * last_success.U(:);

  [S, eig_1, eig_2, eig_3, Eig_1, Eig_2, Eig_3, EIG_1, EIG_2, EIG_3, ...
   sigma_1, sigma_2, sigma_3] = constitutive_problem(E_final, last_success.Ep_prev);
  Sderiv = stiffness_matrix(eig_1, eig_2, eig_3, Eig_1, Eig_2, Eig_3, ...
                            EIG_1, EIG_2, EIG_3, sigma_1, sigma_2, sigma_3);

  return_type = classify_return_types(eig_1, eig_2, eig_3);

  phi = asin(sin_phi);
  cos_phi = sqrt(max(0.0, 1.0 - sin_phi * sin_phi));
  cohesion = c_bar / (2.0 * cos_phi);
  poisson = (3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear));
  young = 2.0 * shear * (1.0 + poisson);

  metadata = struct();
  metadata.nt = nt;
  metadata.nv = nv;
  metadata.accepted_steps = last_success.accepted_steps;
  metadata.final_printed_step = last_success.printed_step;
  metadata.final_zeta = last_success.zeta;
  metadata.final_settlement = last_success.settlement;
  metadata.last_iter = last_success.iter;
  metadata.final_d_zeta = last_success.d_zeta_after_accept;
  metadata.exit_reason_code = exit_reason_code;
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
  metadata.h = h;
  metadata.tolerance = tolerance;
  metadata.it_max = it_max;
  metadata.step_max = step_max;
  metadata.d_zeta_min = d_zeta_min;
  metadata.d_zeta_init = d_zeta_init;
  metadata.d_settle = d_settle;
  metadata.settle_max = settle_max;

  write_metadata(fullfile(outdir, 'meta.txt'), metadata);

  write_ascii_matrix(fullfile(matlab_dir, 'zeta_hist.txt'), zeta_hist);
  write_ascii_matrix(fullfile(matlab_dir, 'settle_hist.txt'), settle_hist);
  write_ascii_matrix(fullfile(matlab_dir, 'U_final.txt'), last_success.U);
  write_ascii_matrix(fullfile(matlab_dir, 'Ep_prev.txt'), last_success.Ep_prev);
  write_ascii_matrix(fullfile(matlab_dir, 'Ep_final_matlab.txt'), last_success.Ep_final);
  write_ascii_matrix(fullfile(matlab_dir, 'E_final.txt'), E_final);
  write_ascii_matrix(fullfile(matlab_dir, 'S_matlab.txt'), S);
  write_ascii_matrix(fullfile(matlab_dir, 'eig_matlab.txt'), [eig_1; eig_2; eig_3]);
  write_ascii_matrix(fullfile(matlab_dir, 'sigma_matlab.txt'), [sigma_1; sigma_2; sigma_3]);
  write_ascii_matrix(fullfile(matlab_dir, 'return_type_matlab.txt'), return_type);
  write_ascii_matrix(fullfile(matlab_dir, 'proj1_matlab.txt'), Eig_1);
  write_ascii_matrix(fullfile(matlab_dir, 'proj2_matlab.txt'), Eig_2);
  write_ascii_matrix(fullfile(matlab_dir, 'proj3_matlab.txt'), Eig_3);
  write_ascii_matrix(fullfile(matlab_dir, 'hess1_matlab_reduced.txt'), EIG_1);
  write_ascii_matrix(fullfile(matlab_dir, 'hess2_matlab_reduced.txt'), EIG_2);
  write_ascii_matrix(fullfile(matlab_dir, 'hess3_matlab_reduced.txt'), EIG_3);
  write_ascii_matrix(fullfile(matlab_dir, 'Sderiv_matlab_reduced.txt'), Sderiv);

  save(fullfile(outdir, 'snapshot.mat'), ...
       'metadata', 'zeta_hist', 'settle_hist', 'E_final', 'S', ...
       'eig_1', 'eig_2', 'eig_3', 'Eig_1', 'Eig_2', 'Eig_3', ...
       'EIG_1', 'EIG_2', 'EIG_3', 'sigma_1', 'sigma_2', 'sigma_3', ...
       'Sderiv', 'return_type');

end

function return_type = classify_return_types(eig_1, eig_2, eig_3)

  global nt
  global c_bar sin_phi lame shear

  trace_E = eig_1 + eig_2 + eig_3;
  f_tr = 2 * shear * ((1 + sin_phi) * eig_1 - (1 - sin_phi) * eig_3) + ...
         2 * lame * sin_phi * trace_E - c_bar;
  gamma_sl = (eig_1 - eig_2) / (1 + sin_phi);
  gamma_sr = (eig_2 - eig_3) / (1 - sin_phi);
  gamma_la = (eig_1 + eig_2 - 2 * eig_3) / (3 - sin_phi);
  gamma_ra = (2 * eig_1 - eig_2 - eig_3) / (3 + sin_phi);

  denom_s = 4 * lame * sin_phi^2 + ...
            2 * shear * (1 + sin_phi)^2 + ...
            2 * shear * (1 - sin_phi)^2;
  denom_l = 4 * lame * sin_phi^2 + ...
            shear * (1 + sin_phi)^2 + ...
            2 * shear * (1 - sin_phi)^2;
  denom_r = 4 * lame * sin_phi^2 + ...
            2 * shear * (1 + sin_phi)^2 + ...
            shear * (1 - sin_phi)^2;

  lambda_s = f_tr / denom_s;
  lambda_l = (shear * ((1 + sin_phi) * (eig_1 + eig_2) - 2 * (1 - sin_phi) * eig_3) + ...
              2 * lame * sin_phi * trace_E - c_bar) / denom_l;
  lambda_r = (shear * (2 * (1 + sin_phi) * eig_1 - (1 - sin_phi) * (eig_2 + eig_3)) + ...
              2 * lame * sin_phi * trace_E - c_bar) / denom_r;

  return_type = 4 * ones(1, nt);

  test_el = (f_tr <= 0);
  return_type(test_el) = 0;

  test_s = (lambda_s <= min(gamma_sl, gamma_sr)) & (~test_el);
  return_type(test_s) = 1;

  test_l = (gamma_sl < gamma_sr) & (lambda_l >= gamma_sl) & ...
           (lambda_l <= gamma_la) & (~(test_el | test_s));
  return_type(test_l) = 2;

  test_r = (gamma_sl > gamma_sr) & (lambda_r >= gamma_sr) & ...
           (lambda_r <= gamma_ra) & (~(test_el | test_s | test_l));
  return_type(test_r) = 3;

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

function ensure_directory(path)

  if exist(path, 'dir') ~= 7
      mkdir(path);
  end

end
