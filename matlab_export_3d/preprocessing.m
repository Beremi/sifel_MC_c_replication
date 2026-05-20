function preprocessing

  global young poisson shear bulk specific_weight
  global x1 x2 x3 y1 y2 z N_h
  global nv nt nint nstrain
  global coord elem surf Dir
  global B WEIGHT F1 K_elast iD jD vD_weight

  [coord, elem, surf, Dir] = regular_mesh_3d(N_h, x1, x2, x3, y1, y2, z);

  nv = size(coord, 2);
  nt = size(elem, 2);
  nint = nt;
  nstrain = 6;

  [B, WEIGHT] = strain_displacement_p1(elem, coord);

  shear_ip = shear * ones(1, nint);
  bulk_ip = bulk * ones(1, nint);
  ELAST = elastic_operator_3d(shear_ip, bulk_ip);

  AUX = reshape(1:nstrain*nint, nstrain, nint);
  iD = repmat(AUX, nstrain, 1);
  jD = kron(AUX, ones(nstrain, 1));
  vD_weight = repmat(WEIGHT, nstrain^2, 1);

  vD_elast = vD_weight .* ELAST;
  D_elast = sparse(iD(:), jD(:), vD_elast(:), nstrain*nint, nstrain*nint);
  K_elast = B' * D_elast * B;
  K_elast = (K_elast + K_elast') / 2;

  f_V_int = [zeros(1, nint); -specific_weight * ones(1, nint); zeros(1, nint)];
  F1 = vector_volume_p1(elem, f_V_int, WEIGHT);

end

function [B, WEIGHT] = strain_displacement_p1(elem, coord)

  n_n = size(coord, 2);
  n_e = size(elem, 2);
  n_p = 4;
  nstrain = 6;

  iB = zeros(nstrain * 3 * n_p, n_e);
  jB = zeros(nstrain * 3 * n_p, n_e);
  vB = zeros(nstrain * 3 * n_p, n_e);
  WEIGHT = zeros(1, n_e);

  for e = 1:n_e
    nodes = elem(:, e);
    P = coord(:, nodes);
    J = [P(:,2) - P(:,1), P(:,3) - P(:,1), P(:,4) - P(:,1)];
    detJ = det(J);
    volume = abs(detJ) / 6;
    invJ = inv(J);

    grad_hat = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];
    grad = grad_hat * invJ;
    WEIGHT(e) = volume;

    row_base = nstrain * (e - 1);
    k = 0;
    for a = 1:n_p
      node = nodes(a);
      dx = grad(a, 1);
      dy = grad(a, 2);
      dz = grad(a, 3);

      rows = row_base + (1:nstrain)';
      cols = [3*node - 2; 3*node - 1; 3*node];
      local = [
        dx, 0,  0
        0,  dy, 0
        0,  0,  dz
        dy, dx, 0
        dz, 0,  dx
        0,  dz, dy
      ];

      for c = 1:3
        idx = k + (1:nstrain);
        iB(idx, e) = rows;
        jB(idx, e) = cols(c);
        vB(idx, e) = local(:, c);
        k = k + nstrain;
      end
    end
  end

  B = sparse(iB(:), jB(:), vB(:), nstrain*n_e, 3*n_n);

end

function ELAST = elastic_operator_3d(shear, bulk)

  iota = [1; 1; 1; 0; 0; 0];
  VOL = iota * iota';
  DEV = diag([1, 1, 1, 1/2, 1/2, 1/2]) - VOL / 3;
  ELAST = 2 * DEV(:) * shear + VOL(:) * bulk;

end

function F = vector_volume_p1(elem, f_V_int, WEIGHT)

  n_n = max(elem(:));
  n_e = size(elem, 2);
  F = zeros(3, n_n);

  for e = 1:n_e
    nodes = elem(:, e);
    for a = 1:4
      F(:, nodes(a)) = F(:, nodes(a)) + WEIGHT(e) * f_V_int(:, e) / 4;
    end
  end

end
