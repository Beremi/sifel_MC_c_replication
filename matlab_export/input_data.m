% ************************************************************************
%
%   input_data  -  input data to the program
%
% ************************************************************************
%

function  input_data
  
  global shear bulk lame       % elastic material parameters
  global c_bar sin_phi         % inelastic material parameters
  global specific_weight       % specific weight
  global x1 x2 x3 y1 y2 h      % geometric and discretization parameters
  global tolerance it_max      % parameters of the Newton method
  global step_max d_zeta_min d_zeta_init d_settle % load process parameters
  global settle_max            % maximal settlement at A
  global DEV Dev VOL iota ELAST
  global Elast Inv_ELAST IDENT Ident % constitutive operators.

% ======================================================================

%
% input data
%

% elastic material parameters
  young = 20000 ;                 % Young's modulus
  poisson = 0.49 ;                 % Poisson's ratio
  
% inelastic material parameters 
  c = 50 ;                         % cohesion
  phi = pi/9;                      % frictional angle  
  
% external forces - gravity load
  specific_weight = 20 ;           % specific weight
 
% geometry of the domain
  x1 = 35 ;                        % length of domain in front of the slope
  x2 = 10 ;                        % length of the slope
  x3 = 30 ;                        % length of domain behind the slope
  y1 = 30 ;                        % height of domain under the slope
  y2 = 10 ;                        % height of the slope
  
% discretization parameter 
  h = 1 ;                          % length of triangles

% parameters of the Newton method
  tolerance = 1e-12;               % tolerance for the Newton method
  it_max = 50;                     % maximal number of Newton iterations
   
% parameters of the loading process  
  step_max = 2000;                 % maximal number of load steps
  d_zeta_min = 1e-5;               % minimal load increment
  d_zeta_init = 1e-1;              % initial load increment
  d_settle = 5e-1;                 % increment of the control variable
  settle_max = 3;                  % maximal settlement at the point A
  
%
% Computed auxilliary parameters
%

% parameters of the Mohr-Coulomb model
  shear   = young/((1+poisson)*2)  ; % shear modulus
  bulk    = young/(3*(1-2*poisson)); % bulk modulus
  lame    = bulk-2*shear/3;          % lame's coefficient (lambda)
  c_bar   = 2*c*cos(phi);             
  sin_phi = sin(phi);

% linear constitutive operators
  IDENT = [1 0  0  0
           0 1  0  0
           0 0 1/2 0
           0 0  0  1];
  DEV = [ 2/3  -1/3   0   -1/3
         -1/3   2/3   0   -1/3
           0     0   1/2    0
         -1/3  -1/3   0    2/3] ;  % deviatoric operator
  Ident = IDENT(1:3,1:3);   
  Dev = DEV(1:3,1:3);              % reduced deviatoric operator
  iota = [1;1;0;1] ;               % unit tensor
  VOL = iota*iota' ;               % volumetric operator
  ELAST = 2*shear*DEV+bulk*VOL;    % elastic operator
  Elast = ELAST(1:3,1:3);          % reduced elastic operator
  Inv_ELAST = inv(ELAST);          % inverse of the elastic operator
  
end 
  


