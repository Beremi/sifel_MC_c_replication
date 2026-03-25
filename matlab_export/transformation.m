% ************************************************************************
%
%  transformation
%
%      input : Q_elem (1*pt) constant values of a quantity Q on elements
%      output: Q_node (1*pu) avarage values of the quantity Q at vertices
%
% ************************************************************************
%
%  Copyright (C) 2015  M. Cermak
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
%*************************************************************************
%

function Q_node=transformation(Q_elem)

 global nt nv area elem
%------------------------------------------------------------------------- 
%
 Q_node=zeros(1,nv);    %
 Area_node=zeros(1,nv); % area of neighbouring elements around a vertex
 for j=1:nt
    Area_node(elem(:,j))=Area_node(elem(:,j))+area(j);
    Q_node(elem(:,j))=Q_node(elem(:,j))+area(j)*Q_elem(j);
 end
 Q_node=Q_node./Area_node;

end