% ************************************************************************
%
%  description of the experimental code 
%
% ************************************************************************
%
%  Copyright (C) 2015  S. Sysala
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% ************************************************************************
%  
%  description of the problem: 
%    - slope stability experiment
%    - plain strain problem
%    - geometry and material parameters from the book 
%       "de Souza Neto, Peric, Owen: Computational Plasticity ..."
%    - perfect plasticity at small strain
%    - Mohr-Coulomb yield criterion
%    - associative plastic flow rule
%    - loading process controled by the gravity load factor, \zeta
%    - loading path: dependence between the settlement of the slope, \alpha
%                    and the gravity load factor, \zeta
%    - implicit Euler time discretization 
%    - return-mapping scheme
%    - linear simplician elements (triangles)
%    - improved return-mapping scheme suggested in 
%        "S. Sysala, M. Cermak, T. Koudelka, J. Kruis, J. Zeman, 
%         R. Blaheta: An improved return-mapping scheme for nonsmooth 
%         plastic potentials: PART II - the Mohr Coulomb yield function.
%         arXive 1508.07435"
%    - semismooth Newton method (nonlinear solver)
%    - vectorization suggested by P. Byczanski
%
%  main program: the file "loading_process.m"
%    - the command "loading_process" triggers the computation
%
%  other functions:
%    - input_data: all input data are introduced in "input_data.m"
%    - regular_mesh: to generate a regular mesh, Dirichlet conditions
%    - preprocessing: to generate strain-displacement relation, areas of 
%                     elements, load vector, elastic stiffness matrix, etc.
%    - newton: the semismooth Newton solver
%    - consitutive problem: a solution of the constitutive problem
%    - transformation: it transforms values of a quantity Q on elements
%                      into avarage values of the quantity Q at vertices
%                      (for purposes of postprocessing).
%
%  global variables:
%    shear bulk            % elastic material parameters
%    c phi                 % inelastic material parameters
%    specific_weight       % specific weight
%    x1 x2 x3 y1 y2 h      % geometric and discretization parameters
%    tolerance it_max alpha_max      % parameters of the Newton method
%    step_max d_zeta_min d_zeta_init d_alpha % load process parameters 
%    DEV Dev VOL iota ELAST Elast Inv_ELAST  % constitutive operators 
%    nt nv        % numbers of triangles and vertices
%    coord        %   2*nv   coordinates of vertices
%    elem         %   3*nt   vertices that create elements 
%    Dir          %   2*nv   logical array for Dirichlet conditions
%    area         %   1*nt   areas of elements
%    EU           % 3nt*2nv  strain-displacement operator
%    F1           %   2*nv   unit load vector 
%    K_elast      % 2nv*2nv  elastic stiffness matrix
%    iMALI jMALI  % auxilliary arrays for stiffness operators
%
%.
function  readme
