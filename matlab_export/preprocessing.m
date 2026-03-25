% ************************************************************************
%
%   preprocessing  - creates a mesh, load vector, elastic stiffness matrix,
%                    strain-displacement matrix, etc.
%
% ************************************************************************
%

function  preprocessing
 
  global specific_weight    % specific weight
  global x1 x2 x3 y1 y2     % geometry of the domain
  global h                  % discretization parameter
  global Elast              % reduced elastic operator
 
  global  nv nt  % number of vertices and triangles
  global  coord  % 2*nv    coordinates of nodes
  global  elem   % 3*nt    elements 
  global  Dir    % 2*nv  logical array - Dirichlet conditions

  global  area         %   1*nt   areas of triangles
  global  EU           % 3nt*2nv  strain-displacement operator
  global  F1           %   2*nv   unit load vector 
  global  K_elast      % 2nv*2nv  elastic stiffness matrix
  global  iMALI jMALI  % auxilliary arrays for stiffness operators

% ======================================================================
%.

% regular triangular mesh 
  regular_mesh         % outputs: coord, elem, nv, nt, Dir

% areas of triangles (area)
% strain-displacement operator (EU)
% unit load vector (F1)
% elastic stiffness matrix (K, iMALI, jMALI)

  area = zeros(1,nt) ;   % areas of triangles
  iEU  = zeros(6,3,nt) ; % indices "i" for the strain-displacement matrix
  jEU  = zeros(6,3,nt) ; % indices "i" for the strain-displacement matrix
  hEU  = zeros(6,3,nt) ; % values of the strain-displacement matrix
  F1   = zeros(2,nv) ;   % unit load vector
  iMALI = zeros(9,nt) ;  % indices "i" for the elastic stiffness matrix
  jMALI = zeros(9,nt) ;  % indices "j" for the elastic stiffness matrix
  hMALI = zeros(9,nt) ;  % values of the elastic stiffness matrix
 
  for  jt = 1:nt      % for any triangle
    
    % area and strain-displacement on the triangle 
    ABC = elem(:,jt) ; % nodal indices of the investigeted triangle
    P   = coord(:,ABC) ; % their coordinates
    MAT = [ P ; ones(1,3) ] ;
    area(jt) = det(MAT)/2 ;   % area of jt-th element (triangle)
          IMA  = inv(MAT)   ;
    DER = IMA(:,1:2);         % derivatives of basis functions
    ist = 3*(jt-1) + (1:3)' ; % auxilliary indices
    for  k = 1:3  % for any node of element
      Vertex = ABC(k) ;% a particular node V
      eu = [ DER(k,1) 0
             0        DER(k,2) 
             DER(k,2) DER(k,1)] ;% strain-displacement operator for V
      iEU(:,k,jt) = [               ist ; ist           ] ;
      jEU(:,k,jt) = [ (2*Vertex-1)*ones(3,1) ; 2*Vertex*ones(3,1) ] ;
      hEU(:,k,jt) = eu(:)                                 ;
    end
    
    % gravity (volume) forces within the element
    F1(2, ABC) = F1(2, ABC)-area(jt)/3 ;     
    
    % elastic stiffness matrix within the element
    iMALI(:,jt) = 3*(jt-1) + [ 1 2 3 1 2 3 1 2 3 ]' ;
    jMALI(:,jt) = 3*(jt-1) + [ 1 1 1 2 2 2 3 3 3 ]' ;
    hMALI(:,jt) = area(jt)*Elast(:)                ;
    
  end % end of for-cycle through the elements

  % strain-displacement operator
  EU = sparse( iEU(:),jEU(:),hEU(:) , 3*nt,2*nv ) ;
  
  % unit load vector
  F1 = specific_weight*F1; 

  % elastic stiffness matrix  
  MALI = sparse( iMALI(:),jMALI(:),hMALI(:) , 3*nt,3*nt ) ;
  K_elast = EU'*MALI*EU ; 
  K_elast = (K_elast'+K_elast)/2 ;     
  
end  % end of function "preprocessing"


  
