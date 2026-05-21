% ************************************************************************
%
%   preprocessing  - creates a mesh, load vector, elastic stiffness matrix,
%                    strain-displacement matrix, etc.
%
% ************************************************************************
%

function  preprocessing

  global specific_weight    % specific weight
  global x1 x2 x3 y1 y2 z   % geometry of the domain
  global N_h                % discretization parameter
  global Elast              % elastic operator

  global  nv nt nint nstrain % number of vertices and tetrahedra
  global  coord              % 3*nv coordinates of nodes
  global  elem               % 4*nt elements
  global  surf               % 3*ns boundary triangles
  global  Dir                % 3*nv logical array - Dirichlet conditions

  global  volume       %   1*nt   volumes of tetrahedra
  global  WEIGHT       %   1*nt   integration weights
  global  EU           % 6nt*3nv  strain-displacement operator
  global  F1           %   3*nv   unit load vector
  global  K_elast      % 3nv*3nv  elastic stiffness matrix
  global  iMALI jMALI  % auxilliary arrays for stiffness operators

% ======================================================================
%.

% regular tetrahedral mesh
  [coord, elem, surf, Dir] = regular_mesh(N_h,x1,x2,x3,y1,y2,z);

  nv = size(coord,2);
  nt = size(elem,2);
  nint = nt;
  nstrain = 6;

% volumes of tetrahedra (volume)
% strain-displacement operator (EU)
% unit load vector (F1)
% elastic stiffness matrix (K, iMALI, jMALI)

  volume = zeros(1,nt);    % volumes of tetrahedra
  iEU  = zeros(18,4,nt);   % indices "i" for the strain-displacement matrix
  jEU  = zeros(18,4,nt);   % indices "j" for the strain-displacement matrix
  hEU  = zeros(18,4,nt);   % values of the strain-displacement matrix
  F1   = zeros(3,nv);      % unit load vector
  iMALI = zeros(36,nt);    % indices "i" for the stiffness operators
  jMALI = zeros(36,nt);    % indices "j" for the stiffness operators
  hMALI = zeros(36,nt);    % values of the elastic stiffness matrix

  mrow = repmat((1:6)',6,1);
  mcol = kron((1:6)',ones(6,1));

  for  jt = 1:nt      % for any tetrahedron

    % volume and strain-displacement on the tetrahedron
    ABC = elem(:,jt);      % nodal indices of the investigated tetrahedron
    P   = coord(:,ABC);    % their coordinates
    MAT = [ P ; ones(1,4) ];
    volume(jt) = abs(det(MAT))/6;
    IMA  = inv(MAT);
    DER = IMA(:,1:3);      % derivatives of basis functions
    ist = 6*(jt-1) + (1:6)';

    for  k = 1:4  % for any node of element
      Vertex = ABC(k);     % a particular node V
      eu = [ DER(k,1) 0        0
             0        DER(k,2) 0
             0        0        DER(k,3)
             DER(k,2) DER(k,1) 0
             0        DER(k,3) DER(k,2)
             DER(k,3) 0        DER(k,1) ]; % upstream strain ordering
      iEU(:,k,jt) = [ ist ; ist ; ist ];
      jEU(:,k,jt) = [ (3*Vertex-2)*ones(6,1)
                      (3*Vertex-1)*ones(6,1)
                       3*Vertex   *ones(6,1) ];
      hEU(:,k,jt) = eu(:);
    end

    % gravity (volume) forces within the element
    F1(2, ABC) = F1(2, ABC)-volume(jt)/4;

    % elastic stiffness matrix within the element
    iMALI(:,jt) = 6*(jt-1) + mrow;
    jMALI(:,jt) = 6*(jt-1) + mcol;
    hMALI(:,jt) = volume(jt)*Elast(:);

  end % end of for-cycle through the elements

  % strain-displacement operator
  EU = sparse( iEU(:),jEU(:),hEU(:) , 6*nt,3*nv );

  % unit load vector
  F1 = specific_weight*F1;

  % elastic stiffness matrix
  MALI = sparse( iMALI(:),jMALI(:),hMALI(:) , 6*nt,6*nt );
  K_elast = EU'*MALI*EU;
  K_elast = (K_elast'+K_elast)/2;

  WEIGHT = volume;

end
