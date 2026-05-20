function [U, Ep, iter1, last_response] = newton(U_init, Ep_init, F_ext)

  global young poisson cohesion phi
  global B WEIGHT Dir K_elast iD jD vD_weight
  global nint nstrain
  global tolerance it_max settle_max

  U = U_init;
  dU = zeros(size(U));
  iter1 = 0;

  [E_sifel, S_sifel, DS_sifel, return_type, Ep, eig_values, sigma_values] = ...
    evaluate_material(U, Ep_init);

  while true
    iter1 = iter1 + 1;

    if iter1 >= it_max
      warning('Maximal number of Newton iterations was achieved.');
      break
    end

    S_internal = sifel_to_internal_rows(S_sifel);
    DS_internal = sifel_to_internal_tangent(DS_sifel);

    F_int = B' * reshape(repmat(WEIGHT, nstrain, 1) .* S_internal, [], 1);
    F_int = reshape(F_int, 3, []);

    vD = vD_weight .* DS_internal;
    D_p = sparse(iD(:), jD(:), vD(:), nstrain*nint, nstrain*nint);
    K_tangent = B' * D_p * B;
    K_tangent = (K_tangent + K_tangent') / 2;

    residual = F_ext - F_int;
    dU(:) = 0;
    dU(Dir) = K_tangent(Dir, Dir) \ residual(Dir);

    U_new = U + dU;

    q1 = sqrt(abs(dU(:)' * K_elast * dU(:)));
    q2 = sqrt(abs(U(:)' * K_elast * U(:)));
    q3 = sqrt(abs(U_new(:)' * K_elast * U_new(:)));
    criterion = q1 / (q2 + q3 + eps);

    U = U_new;
    [E_sifel, S_sifel, DS_sifel, return_type, Ep, eig_values, sigma_values] = ...
      evaluate_material(U, Ep_init);

    if criterion < tolerance
      break
    end

    if max(abs(U(:))) > 2 * settle_max
      warning('Too large displacement in the Newton solver.');
      iter1 = it_max;
      break
    end
  end

  last_response = struct();
  last_response.E_sifel = E_sifel;
  last_response.S_sifel = S_sifel;
  last_response.DS_sifel = DS_sifel;
  last_response.return_type = return_type;
  last_response.eig_values = eig_values;
  last_response.sigma_values = sigma_values;

end

function [E_sifel, S_sifel, DS_sifel, return_type, Ep, eig_values, sigma_values] = evaluate_material(U, Ep_prev)

  global young poisson cohesion phi B nint nstrain

  perm = [1, 2, 3, 6, 5, 4];
  E_internal = reshape(B * U(:), nstrain, nint);
  E_sifel = E_internal(perm, :);

  [S_sifel, DS_sifel, return_type, Ep, eig_values, sigma_values] = ...
    constitutive_problem_3d(E_sifel, Ep_prev, young, poisson, cohesion, phi);

end

function A_internal = sifel_to_internal_rows(A_sifel)
  perm = [1, 2, 3, 6, 5, 4];
  A_internal = A_sifel(perm, :);
end

function DS_internal = sifel_to_internal_tangent(DS_sifel)
  perm = [1, 2, 3, 6, 5, 4];
  nt = size(DS_sifel, 2);
  DS_internal = zeros(size(DS_sifel));

  for ip = 1:nt
    D_sifel = reshape(DS_sifel(:, ip), 6, 6);
    D_internal = D_sifel(perm, perm);
    DS_internal(:, ip) = D_internal(:);
  end
end
