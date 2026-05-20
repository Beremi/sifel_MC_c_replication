function input_data

  global young poisson cohesion phi
  global shear bulk lame c_bar sin_phi
  global specific_weight
  global x1 x2 x3 y1 y2 z N_h
  global tolerance it_max settle_max
  global step_max d_zeta_min d_zeta_init d_settle
  global target_displacement post_plastic_d_zeta
  global N_h_override
  global nstrain

  young = 20000;
  poisson = 0.49;
  cohesion = 50;
  phi = pi / 9;

  specific_weight = 20;

  % Small deterministic 3D homogeneous slope, using the upstream P1 mesh path.
  x1 = 3;
  x2 = 2;
  x3 = 3;
  y1 = 2;
  y2 = 1;
  z = 1;
  N_h = 4;
  if ~isempty(N_h_override)
    N_h = N_h_override;
  end

  tolerance = 1e-10;
  it_max = 50;

  step_max = 500;
  d_zeta_min = 1e-4;
  d_zeta_init = 0.5;
  post_plastic_d_zeta = 5.0;
  target_displacement = 5.0e-1;
  d_settle = 1.0;
  settle_max = 50.0;

  nstrain = 6;

  shear = young / (2 * (1 + poisson));
  bulk = young / (3 * (1 - 2 * poisson));
  lame = bulk - 2 * shear / 3;
  c_bar = 2 * cohesion * cos(phi);
  sin_phi = sin(phi);

end
