% ************************************************************************
%
%   input_data  -  input data to the 3D program
%
% ************************************************************************
%

function  input_data

  global young poisson c cohesion phi
  global shear bulk lame       % elastic material parameters
  global c_bar sin_phi         % inelastic material parameters
  global specific_weight       % specific weight
  global x1 x2 x3 y1 y2 z      % geometric parameters
  global N_h                   % discretization parameter
  global tolerance it_max      % parameters of the Newton method
  global step_max d_zeta_min d_zeta_init d_settle % load process parameters
  global settle_max            % maximal settlement
  global target_displacement post_plastic_d_zeta
  global DEV VOL iota ELAST
  global Elast Inv_ELAST IDENT % constitutive operators
  global N_h_override
  global nstrain

% ======================================================================

%
% input data
%

% elastic material parameters
  young = 20000;                    % Young's modulus
  poisson = 0.49;                   % Poisson's ratio

% inelastic material parameters
  c = 50;                           % cohesion
  cohesion = c;                     % export metadata name
  phi = pi/9;                       % frictional angle

% external forces - gravity load
  specific_weight = 20;             % specific weight

% geometry of the domain
  x1 = 3;                           % length of domain in front of the slope
  x2 = 2;                           % length of the slope
  x3 = 3;                           % length of domain behind the slope
  y1 = 2;                           % height of domain under the slope
  y2 = 1;                           % height of the slope
  z  = 1;                           % width of the 3D slice

% discretization parameter
  N_h = 4;                          % uniform 3D mesh density
  if ~isempty(N_h_override)
    N_h = N_h_override;
  end

% parameters of the Newton method
  tolerance = 1e-10;                % tolerance for the Newton method
  it_max = 50;                      % maximal number of Newton iterations

% parameters of the loading process
  step_max = 500;                   % maximal number of load steps
  d_zeta_min = 1e-4;                % minimal load increment
  d_zeta_init = 5e-1;               % initial load increment
  post_plastic_d_zeta = 5.0;        % load increment after first plasticity
  target_displacement = 5.0e-1;     % target displacement for the export
  d_settle = 1.0;                   % increment of the control variable
  settle_max = 50.0;                % maximal settlement

%
% Computed auxilliary parameters
%

% parameters of the Mohr-Coulomb model
  shear   = young/((1+poisson)*2);  % shear modulus
  bulk    = young/(3*(1-2*poisson));% bulk modulus
  lame    = bulk-2*shear/3;         % Lame's coefficient (lambda)
  c_bar   = 2*cohesion*cos(phi);
  sin_phi = sin(phi);

% linear constitutive operators
  nstrain = 6;
  IDENT = diag([1, 1, 1, 1/2, 1/2, 1/2]);
  iota = [1; 1; 1; 0; 0; 0];       % unit tensor
  VOL = iota*iota';                % volumetric operator
  DEV = IDENT - VOL/3;             % deviatoric operator
  ELAST = 2*shear*DEV + bulk*VOL;  % elastic operator
  Elast = ELAST;
  Inv_ELAST = inv(ELAST);          % inverse elastic operator

end
