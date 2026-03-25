% ************************************************************************
%
%  regular_mesh  -  create a regular mesh
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
%                      x1               
%              -------------------A
%             |                   \      
% zeroth      |                    \ x2    
% displacement| y2                  \    
% at x-       |                      \   
% direction   |                       \      x3 
%             |                        ------------|
%             |                                    | zeroth displacement
%             |                                 y1 | at x-direction
%             |                                    |
%             |                                    | 
%             --------------------------------------
%             zeroth displacement at both direction
%                     
% ************************************************************************
%

function  regular_mesh
 
  global  x1 x2 x3 y1 y2   % geometry of the domain
  global  h                % diskretization parameter

  global  coord %  2*pu  coordinates 
  global  elem  %  3*nt  elements
  global  nv nt %  number of nodes and elements
  global  Dir   %  2*pu logical array that detects the Dirichlet conditions

% ======================================================================

% number of intervals in x direction
 nx1 = round(x1/h) ;
 nx2 = round(x2/h) ;
 nx3 = round(x3/h) ;
 nx = nx1 + nx2 + nx3 ;

% number of intervals in y direction
 ny1 = round(y1/h) ;
 ny2 = nx2 ;
 ny = ny1 + ny2 ;

% lengths of intervals in x and y directions
 hx1 = x1/nx1 ;
 hx2 = x2/nx2 ;
 hx3 = x3/nx3 ;
 hy1 = y1/ny1 ;
 hy2 = hx2*y2/x2 ;
 
% coordinates in x direction
 coord_x = zeros(1,nx+1) ;
 for i=1:nx1
     coord_x(i+1) = coord_x(i) + hx1 ;
 end
 for i=(nx1+1):(nx1+nx2)
     coord_x(i+1) = coord_x(i) + hx2 ;
 end
 for i=(nx1+nx2+1):nx
     coord_x(i+1) = coord_x(i) + hx3 ;
 end
 
% coordinates in y direction
 coord_y = zeros(1,ny+1) ;
 for j=1:ny1
     coord_y(j+1) = coord_y(j) + hy1 ;
 end
 for j=(ny1+1):ny
     coord_y(j+1) = coord_y(j) + hy2 ;
 end

% indexing of nodes and their coordinates
 V = zeros(nx+1,ny+1) ;     % indices of vertices in a matrix ordering
 coord = zeros(2,(nx+1)*(ny+1)) ;
 nv = 0 ;
 for j=1:(ny1+1)
     for i=1:(nx+1)
         nv = nv+1 ;
         V(i,j) = nv ;
         coord(:,nv) = [coord_x(i); coord_y(j)] ;
     end
 end
 for j=1:ny2
     for i=1:(nx1+nx2+1-j)
         nv = nv+1 ;
         V(i,ny1+1+j) = nv ;
         coord(:,nv) = [coord_x(i); coord_y(ny1+1+j)] ;
     end
 end
 coord = coord(:,1:nv);
 
% elements

 % first alternative - diagonals from left-up to right-down
 elem = zeros(3,2*nx*ny) ;
 nt = 0 ;
 for  j = 1:ny1
   for  i = 1:nx
      nt = nt + 1 ;
      elem(:,nt) = [ V(i,j)   ; V(i+1,j)   ; V(i,j+1) ] ;
      nt = nt + 1 ;
      elem(:,nt) = [ V(i,j+1) ; V(i+1,j)   ; V(i+1,j+1) ] ;
   end
 end
 for  j = 1:ny2
   for  i = 1:(nx1+nx2-j)
      nt = nt + 1 ;
      elem(:,nt) = [ V(i,ny1+j)   ; V(i+1,ny1+j)   ; V(i,ny1+j+1) ] ;
      nt = nt + 1 ;
      elem(:,nt) = [ V(i,ny1+j+1) ; V(i+1,ny1+j)   ; V(i+1,ny1+j+1) ] ;
   end
   i = i+1 ;
   nt = nt + 1 ;
   elem(:,nt) = [ V(i,ny1+j)   ; V(i+1,ny1+j)   ; V(i,ny1+j+1) ] ;
 end
 elem = elem(:,1:nt) ;

%  % second alternative - diagonals from left-down to right-up 
%  elem = zeros(3,2*nx*ny) ;
%  nt = 0 ;
%  for  j = 1:ny1
%    for  i = 1:nx
%       nt = nt + 1 ;
%       elem(:,nt) = [ V(i,j)   ; V(i+1,j)   ; V(i+1,j+1) ] ;
%       nt = nt + 1 ;
%       elem(:,nt) = [ V(i+1,j+1) ; V(i,j+1)   ; V(i,j) ] ;
%    end
%  end
%  for  j = 1:ny2
%    for  i = 1:(nx1+nx2-j)
%       nt = nt + 1 ;
%       elem(:,nt) = [ V(i,ny1+j)   ; V(i+1,ny1+j)   ; V(i+1,ny1+j+1) ] ;
%       nt = nt + 1 ;
%       elem(:,nt) = [ V(i+1,ny1+j+1) ; V(i,ny1+j+1)   ; V(i,ny1+j) ] ;
%    end
%    i = i+1 ;
%    nt = nt + 1 ;
%    elem(:,nt) = [ V(i,ny1+j)   ; V(i+1,ny1+j)   ; V(i,ny1+j+1) ] ;
%  end
%  elem = elem(:,1:nt) ;
 
% logical array that detects the Dirichlet conditions
 x_max=max(coord(1,:));
 Dir=false(2,nv);
 Dir(1,:) = (coord(1,:)>0)&(coord(2,:)>0)&(coord(1,:)<x_max) ;
 Dir(2,:) = coord(2,:)>0 ;

 
% % figure with mesh
% figure
%   hold on
%     for  jt = 1:size(elem,2) 
%       P = coord(:,elem(:,jt)) ;
%         e = P(:,3) ;
%       for  j = 1:3
%         f = P(:,j) ;
%           plot( [e(1) f(1)],[e(2) f(2)] )
%         e = f ;
%       end
%     end
%   box on
%   axis equal

end