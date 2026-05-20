function [S_sifel, DS_sifel, return_type, Ep_final_sifel, eig_values, sigma_values] = ...
  constitutive_problem_3d(E_new_sifel, Ep_prev_sifel, young, poisson, cohesion, phi)

  perm = [1, 2, 3, 6, 5, 4];
  nt = size(E_new_sifel, 2);

  shear = young / (2 * (1 + poisson));
  bulk = young / (3 * (1 - 2 * poisson));
  lame = bulk - 2 * shear / 3;
  sin_phi = sin(phi);
  cos_phi = cos(phi);
  c_bar = 2 * cohesion * cos_phi;

  E_new_internal = E_new_sifel(perm, :);
  Ep_prev_internal = Ep_prev_sifel(perm, :);
  E_trial_internal = E_new_internal - Ep_prev_internal;

  S_sifel = zeros(6, nt);
  DS_sifel = zeros(36, nt);
  return_type = zeros(1, nt);
  Ep_final_sifel = zeros(6, nt);
  eig_values = zeros(3, nt);
  sigma_values = zeros(3, nt);

  for ip = 1:nt
    [S_internal, D_internal, ret, eig_i, sigma_i] = ...
      response_internal(E_trial_internal(:, ip), c_bar, sin_phi, shear, bulk, lame);

    S_ip_sifel = S_internal(perm);
    D_ip_sifel = D_internal(perm, perm);
    eps_el_sifel = elastic_strain_from_stress(S_ip_sifel, young, poisson, shear);

    S_sifel(:, ip) = S_ip_sifel;
    DS_sifel(:, ip) = D_ip_sifel(:);
    return_type(ip) = ret;
    Ep_final_sifel(:, ip) = E_new_sifel(:, ip) - eps_el_sifel;
    eig_values(:, ip) = eig_i(:);
    sigma_values(:, ip) = sigma_i(:);
  end

end

function [S, D, ret, eig, sigma] = response_internal(E_trial, c_bar, sin_phi, shear, bulk, lame)

  IDENT_DIAG = [1; 1; 1; 1/2; 1/2; 1/2];
  iota = [1; 1; 1; 0; 0; 0];
  IDENT = diag(IDENT_DIAG);

  E_tr = IDENT * E_trial;
  E_square = [
    E_tr(1)^2 + E_tr(4)^2 + E_tr(6)^2
    E_tr(2)^2 + E_tr(4)^2 + E_tr(5)^2
    E_tr(3)^2 + E_tr(5)^2 + E_tr(6)^2
    E_tr(1)*E_tr(4) + E_tr(2)*E_tr(4) + E_tr(5)*E_tr(6)
    E_tr(4)*E_tr(6) + E_tr(2)*E_tr(5) + E_tr(3)*E_tr(5)
    E_tr(1)*E_tr(6) + E_tr(4)*E_tr(5) + E_tr(3)*E_tr(6)
  ];

  I1 = E_tr(1) + E_tr(2) + E_tr(3);
  I2 = E_tr(1)*E_tr(2) + E_tr(1)*E_tr(3) + E_tr(2)*E_tr(3) - ...
       E_tr(4)^2 - E_tr(5)^2 - E_tr(6)^2;
  I3 = E_tr(1)*E_tr(2)*E_tr(3) - E_tr(3)*E_tr(4)^2 - ...
       E_tr(2)*E_tr(6)^2 - E_tr(1)*E_tr(5)^2 + ...
       2*E_tr(4)*E_tr(5)*E_tr(6);

  Q = max(0, (I1^2 - 3*I2) / 9);
  R = (-2*I1^3 + 9*I1*I2 - 27*I3) / 54;
  theta0 = 0;
  if Q > 0
    theta0 = R / sqrt(Q^3);
  end
  theta = acos(min(max(theta0, -1), 1)) / 3;

  eig = [
    -2*sqrt(Q)*cos(theta + 2*pi/3) + I1/3
    -2*sqrt(Q)*cos(theta - 2*pi/3) + I1/3
    -2*sqrt(Q)*cos(theta) + I1/3
  ];

  f_tr = 2*shear*((1 + sin_phi)*eig(1) - (1 - sin_phi)*eig(3)) + ...
         2*lame*sin_phi*I1 - c_bar;
  gamma_sl = (eig(1) - eig(2)) / (1 + sin_phi);
  gamma_sr = (eig(2) - eig(3)) / (1 - sin_phi);
  gamma_la = (eig(1) + eig(2) - 2*eig(3)) / (3 - sin_phi);
  gamma_ra = (2*eig(1) - eig(2) - eig(3)) / (3 + sin_phi);

  denom_s = 4*lame*sin_phi^2 + 4*shear*(1 + sin_phi^2);
  denom_l = 4*lame*sin_phi^2 + shear*(1 + sin_phi)^2 + 2*shear*(1 - sin_phi)^2;
  denom_r = 4*lame*sin_phi^2 + 2*shear*(1 + sin_phi)^2 + shear*(1 - sin_phi)^2;
  denom_a = 4*bulk*sin_phi^2;

  lambda_s = f_tr / denom_s;
  lambda_l = (shear*((1 + sin_phi)*(eig(1) + eig(2)) - 2*(1 - sin_phi)*eig(3)) + ...
              2*lame*sin_phi*I1 - c_bar) / denom_l;
  lambda_r = (shear*(2*(1 + sin_phi)*eig(1) - (1 - sin_phi)*(eig(2) + eig(3))) + ...
              2*lame*sin_phi*I1 - c_bar) / denom_r;
  lambda_a = (2*bulk*sin_phi*I1 - c_bar) / denom_a; %#ok<NASGU>

  if f_tr <= 0
    ret = 0;
  elseif lambda_s <= min(gamma_sl, gamma_sr)
    ret = 1;
  elseif (gamma_sl < gamma_sr) && (lambda_l >= gamma_sl) && (lambda_l <= gamma_la)
    ret = 2;
  elseif (gamma_sl > gamma_sr) && (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra)
    ret = 3;
  else
    ret = 4;
  end

  S = zeros(6, 1);
  D = zeros(6, 6);
  sigma = zeros(3, 1);
  DER_E_square = der_e_square(E_tr);

  switch ret
    case 0
      sigma = [
        lame*I1 + 2*shear*eig(1)
        lame*I1 + 2*shear*eig(2)
        lame*I1 + 2*shear*eig(3)
      ];
      S = lame*I1*iota + 2*shear*E_tr;
      D = lame*(iota*iota') + 2*shear*IDENT;

    case 1
      Eig_1 = eigenprojection(E_square, E_tr, eig(1), eig(2), eig(3), iota);
      Eig_2 = eigenprojection(E_square, E_tr, eig(2), eig(1), eig(3), iota);
      Eig_3 = eigenprojection(E_square, E_tr, eig(3), eig(1), eig(2), iota);
      sigma = [
        lame*I1 + 2*shear*eig(1) - lambda_s*(2*lame*sin_phi + 2*shear*(1 + sin_phi))
        lame*I1 + 2*shear*eig(2) - lambda_s*(2*lame*sin_phi)
        lame*I1 + 2*shear*eig(3) - lambda_s*(2*lame*sin_phi - 2*shear*(1 - sin_phi))
      ];
      S = sigma(1)*Eig_1 + sigma(2)*Eig_2 + sigma(3)*Eig_3;

      EIG_1 = smooth_eig_derivative(DER_E_square, IDENT, Eig_1, Eig_2, Eig_3, eig(1), eig(2), eig(3));
      EIG_2 = smooth_eig_derivative(DER_E_square, IDENT, Eig_2, Eig_1, Eig_3, eig(2), eig(1), eig(3));
      EIG_3 = smooth_eig_derivative(DER_E_square, IDENT, Eig_3, Eig_1, Eig_2, eig(3), eig(1), eig(2));
      D_phi = 2*shear*((1 + sin_phi)*Eig_1 - (1 - sin_phi)*Eig_3) + 2*lame*sin_phi*iota;
      D = sigma(1)*EIG_1 + sigma(2)*EIG_2 + sigma(3)*EIG_3 + ...
          lame*(iota*iota') + 2*shear*(Eig_1*Eig_1' + Eig_2*Eig_2' + Eig_3*Eig_3') - ...
          (D_phi*D_phi') / denom_s;

    case 2
      Eig_3 = eigenprojection(E_square, E_tr, eig(3), eig(1), eig(2), iota);
      Eig_12 = iota - Eig_3;
      sigma(1) = lame*I1 + shear*(eig(1) + eig(2)) - ...
                 lambda_l*(2*lame*sin_phi + shear*(1 + sin_phi));
      sigma(2) = sigma(1);
      sigma(3) = lame*I1 + 2*shear*eig(3) - ...
                 lambda_l*(2*lame*sin_phi - 2*shear*(1 - sin_phi));
      S = sigma(1)*Eig_12 + sigma(3)*Eig_3;

      EIG_3 = left_eig3_derivative(DER_E_square, IDENT, E_tr, Eig_12, Eig_3, eig(1), eig(2), eig(3));
      D_phi = shear*((1 + sin_phi)*Eig_12 - 2*(1 - sin_phi)*Eig_3) + 2*lame*sin_phi*iota;
      D = (sigma(3) - sigma(1))*EIG_3 + lame*(iota*iota') + ...
          shear*(Eig_12*Eig_12' + 2*(Eig_3*Eig_3')) - (D_phi*D_phi') / denom_l;

    case 3
      Eig_1 = eigenprojection(E_square, E_tr, eig(1), eig(2), eig(3), iota);
      Eig_23 = iota - Eig_1;
      sigma(1) = lame*I1 + 2*shear*eig(1) - ...
                 lambda_r*(2*lame*sin_phi + 2*shear*(1 + sin_phi));
      sigma(3) = lame*I1 + shear*(eig(2) + eig(3)) - ...
                 lambda_r*(2*lame*sin_phi - shear*(1 - sin_phi));
      sigma(2) = sigma(3);
      S = sigma(1)*Eig_1 + sigma(3)*Eig_23;

      EIG_1 = right_eig1_derivative(DER_E_square, IDENT, E_tr, Eig_1, Eig_23, eig(1), eig(2), eig(3));
      D_phi = shear*(2*(1 + sin_phi)*Eig_1 - (1 - sin_phi)*Eig_23) + 2*lame*sin_phi*iota;
      D = (sigma(1) - sigma(3))*EIG_1 + lame*(iota*iota') + ...
          shear*(2*(Eig_1*Eig_1') + Eig_23*Eig_23') - (D_phi*D_phi') / denom_r;

    otherwise
      sigma(:) = c_bar / (2 * sin_phi);
      S = iota * sigma(1);
      D = zeros(6, 6);
  end

end

function eps_el = elastic_strain_from_stress(S, young, poisson, shear)
  eps_el = zeros(6, 1);
  eps_el(1) = (S(1) - poisson*(S(2) + S(3))) / young;
  eps_el(2) = (S(2) - poisson*(S(1) + S(3))) / young;
  eps_el(3) = (S(3) - poisson*(S(1) + S(2))) / young;
  eps_el(4) = S(4) / shear;
  eps_el(5) = S(5) / shear;
  eps_el(6) = S(6) / shear;
end

function Eig = eigenprojection(E_square, E_tr, eig_i, eig_j, eig_k, iota)
  Eig = (E_square - (eig_j + eig_k)*E_tr + iota*(eig_j*eig_k)) / ...
        ((eig_i - eig_j) * (eig_i - eig_k));
end

function DER = der_e_square(E_tr)
  DER = zeros(6, 6);
  DER(:, 1) = [2*E_tr(1); 0; 0; E_tr(4); 0; E_tr(6)];
  DER(:, 2) = [0; 2*E_tr(2); 0; E_tr(4); E_tr(5); 0];
  DER(:, 3) = [0; 0; 2*E_tr(3); 0; E_tr(5); E_tr(6)];
  DER(:, 4) = [E_tr(4); E_tr(4); 0; 0.5*(E_tr(1) + E_tr(2)); 0.5*E_tr(6); 0.5*E_tr(5)];
  DER(:, 5) = [0; E_tr(5); E_tr(5); 0.5*E_tr(6); 0.5*(E_tr(2) + E_tr(3)); 0.5*E_tr(4)];
  DER(:, 6) = [E_tr(6); 0; E_tr(6); 0.5*E_tr(5); 0.5*E_tr(4); 0.5*(E_tr(1) + E_tr(3))];
end

function EIG_i = smooth_eig_derivative(DER, IDENT, Eig_i, Eig_j, Eig_k, eig_i, eig_j, eig_k)
  EIG_i = (DER - IDENT*(eig_j + eig_k) - ...
           (2*eig_i - eig_j - eig_k)*(Eig_i*Eig_i') - ...
           (eig_j - eig_k)*(Eig_j*Eig_j' - Eig_k*Eig_k')) / ...
          ((eig_i - eig_j) * (eig_i - eig_k));
end

function EIG_3 = left_eig3_derivative(DER, IDENT, E_tr, Eig_12, Eig_3, eig_1, eig_2, eig_3)
  EIG_3 = (DER - IDENT*(eig_1 + eig_2) - E_tr*Eig_12' - Eig_12*E_tr' + ...
           (eig_1 + eig_2)*(Eig_12*Eig_12') + ...
           (eig_1 + eig_2 - 2*eig_3)*(Eig_3*Eig_3') + ...
           eig_3*(Eig_12*Eig_3' + Eig_3*Eig_12')) / ...
          ((eig_3 - eig_1) * (eig_3 - eig_2));
end

function EIG_1 = right_eig1_derivative(DER, IDENT, E_tr, Eig_1, Eig_23, eig_1, eig_2, eig_3)
  EIG_1 = (DER - IDENT*(eig_2 + eig_3) - E_tr*Eig_23' - Eig_23*E_tr' + ...
           (eig_2 + eig_3)*(Eig_23*Eig_23') + ...
           (eig_2 + eig_3 - 2*eig_1)*(Eig_1*Eig_1') + ...
           eig_1*(Eig_23*Eig_1' + Eig_1*Eig_23')) / ...
          ((eig_1 - eig_2) * (eig_1 - eig_3));
end
