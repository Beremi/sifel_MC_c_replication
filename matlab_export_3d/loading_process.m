function result = loading_process

  global nv nint
  global F1
  global step_max d_zeta_min d_zeta_init d_settle settle_max
  global target_displacement post_plastic_d_zeta
  global it_max

  input_data;
  preprocessing;

  zeta_hist = zeros(1, step_max);
  settle_hist = zeros(1, step_max);
  max_displacement_hist = zeros(1, step_max);
  iter_hist = zeros(1, step_max);
  plastic_count_hist = zeros(1, step_max);
  return_type_count_hist = zeros(5, step_max);

  d_zeta = d_zeta_init;
  settle_old = 0;
  zeta_old = 0;
  U_old = zeros(3, nv);
  Ep_old = zeros(6, nint);
  step = 1;
  accepted_steps = 0;
  post_plastic_reached = false;
  post_plastic_increment_initialized = false;
  exit_reason_code = 0;

  last_success = struct( ...
    'accepted', false, ...
    'printed_step', 1, ...
    'accepted_steps', 0, ...
    'zeta', 0, ...
    'settlement', 0, ...
    'max_displacement', 0, ...
    'iter', 0, ...
    'd_zeta_after_accept', d_zeta, ...
    'U', zeros(3, nv), ...
    'Ep_prev', zeros(6, nint), ...
    'Ep_final', zeros(6, nint), ...
    'response', struct());

  while true
    zeta = zeta_old + d_zeta;
    F = zeta * F1;
    Ep_prev_step = Ep_old;

    [U, Ep, iter, response] = newton(U_old, Ep_old, F);

    if iter == it_max
      d_zeta = max(d_zeta / 2, d_zeta_min);
    else
      step = step + 1;
      accepted_steps = accepted_steps + 1;
      zeta_old = zeta;
      zeta_hist(step) = zeta;
      settle = max(0, -min(U(2, :)));
      max_displacement = max(sqrt(sum(U.^2, 1)));
      settle_hist(step) = settle;
      max_displacement_hist(step) = max_displacement;
      iter_hist(step) = iter;

      return_counts = histcounts(response.return_type, -0.5:1:4.5);
      plastic_count = sum(response.return_type > 0);
      plastic_count_hist(step) = plastic_count;
      return_type_count_hist(:, step) = return_counts(:);

      if (settle - settle_old) > d_settle
        warning('Too large increment of the settlement.');
        d_zeta = max(d_zeta / 2, d_zeta_min);
      end

      settle_old = settle;
      U_old = U;
      Ep_old = Ep;

      if plastic_count > 0
        post_plastic_reached = true;
        if ~post_plastic_increment_initialized
          d_zeta = post_plastic_d_zeta;
          post_plastic_increment_initialized = true;
        end
      end

      last_success.accepted = true;
      last_success.printed_step = step;
      last_success.accepted_steps = accepted_steps;
      last_success.zeta = zeta;
      last_success.settlement = settle;
      last_success.max_displacement = max_displacement;
      last_success.iter = iter;
      last_success.d_zeta_after_accept = d_zeta;
      last_success.U = U;
      last_success.Ep_prev = Ep_prev_step;
      last_success.Ep_final = Ep;
      last_success.response = response;

      disp([' step=', num2str(step), ', settlement=', ...
            num2str(settle), ', zeta=', num2str(zeta), ', d_zeta=', ...
            num2str(d_zeta), ', iter=', num2str(iter), ...
            ', max_u=', num2str(max_displacement), ...
            ', plastic_points=', num2str(plastic_count)])
    end

    if post_plastic_reached && (last_success.max_displacement >= target_displacement)
      exit_reason_code = 4;
      break
    end

    if settle_old >= settle_max
      warning('Too large settlement in the 3D loading process.');
      exit_reason_code = 1;
      break
    end

    if d_zeta == d_zeta_min
      warning('Too small load increments in the 3D loading process.');
      exit_reason_code = 2;
      break
    end

    if step >= step_max
      warning('Maximal number of load steps was achieved.');
      exit_reason_code = 3;
      break
    end
  end

  if ~last_success.accepted
    error('No successful 3D load step was completed.');
  end

  result = struct();
  result.last_success = last_success;
  result.exit_reason_code = exit_reason_code;
  result.zeta_hist = zeta_hist(1:step);
  result.settle_hist = settle_hist(1:step);
  result.max_displacement_hist = max_displacement_hist(1:step);
  result.iter_hist = iter_hist(1:step);
  result.plastic_count_hist = plastic_count_hist(1:step);
  result.return_type_count_hist = return_type_count_hist(:, 1:step);

end
