% ************************************************************************
%
%  3D MATLAB export
%
% ************************************************************************
%
% The file structure intentionally follows ../matlab_export file by file:
%
%   input_data.m
%   regular_mesh.m
%   preprocessing.m
%   constitutive_problem.m
%   stiffness_matrix.m
%   newton.m
%   loading_process.m
%   transformation.m
%
% The 3D constitutive arrays use the upstream slope_stability ordering
% [xx, yy, zz, xy, yz, xz] inside the MATLAB solver. Export files for the
% C++ SIFEL-shaped comparison are written in SIFEL ordering
% [xx, yy, zz, yz, xz, xy].
