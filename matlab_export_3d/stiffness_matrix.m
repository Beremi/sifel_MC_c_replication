% ************************************************************************
%
%  stiffness_matrix   -  assembling of the 3D tangent stiffness matrix for
%                        the Mohr-Coulomb problem
%
% ************************************************************************

function   Sderiv          ... % 36*nt consistent tangent stiffness matrix
             =             ...
  stiffness_matrix...
     (DS)                      % tangent from the 3D constitutive problem

% .........................................................................

  Sderiv = DS;

end
