% ************************************************************************
%
%  loading_process - the main program, for more details: readme.m
%
% ************************************************************************
%
  
% global variables defined in "input_data"
  global shear bulk lame       % elastic material parameters
  global c_bar sin_phi         % inelastic material parameters
  global specific_weight       % specific weight
  global x1 x2 x3 y1 y2 h      % geometric and discretization parameters
  global tolerance it_max settle_max     % parameters of the Newton method
  global step_max d_zeta_min d_zeta_init d_settle % load process parameters 
  global DEV Dev VOL iota ELAST 
  global Elast Inv_ELAST IDENT Ident % constitutive operators
  
% global variables defined in "preprocessing" 
  global nt nv        % numbers of triangles and vertices
  global coord        %   2*nv   coordinates of vertices
  global elem         %   3*nt   vertices that create elements 
  global Dir          %   2*nv   logical array for Dirichlet conditions
  global area         %   1*nt   areas of elements
  global EU           % 3nt*2nv  strain-displacement operator
  global F1           %   2*nv   unit load vector 
  global K_elast      % 2nv*2nv  elastic stiffness matrix
  global iMALI jMALI  % auxilliary arrays for stiffness operators

   
% ========================================================================

%
% input data
%
  input_data     % outputs: shear,bulk,lame,c_bar,sin_phi, specific_weight,
                 %          x1, x2, x3, y1, y2, h, tolerance, it_max,
                 %          step_max, d_zeta_min, d_zeta_init, d_alpha,
                 %          DEV, VOL, ELAST, iota, Dev, Elast, inv_ELAST
                 %          IDENT, Ident

%
% preprocessing 
%
  preprocessing  % outputs: coord, elem, area, EU, F1, Dir, nt, nv, 
                 %          iMALI, jMALI, K_elast 

%
% initialization 
%
  zeta_hist=zeros(1,step_max);  % history of the gravity load factor
  settle_hist=zeros(1,step_max); % history of the settlement at A
  d_zeta=d_zeta_init;           % increment of the gravity load factor
  settle_old=0;
  zeta_old=0;
  U_old=zeros(2,nv);            % zeroth initial conditions
  Ep_old=zeros(4,nt);           % zeroth initial conditions
  step=1;

%
% loading process
%
  while true   
      
      zeta=zeta_old+d_zeta;  % investigated gravity load factor
      F=zeta*F1;             % the corresponding load vector
      
      % damped Newton solver
      [U, Ep, iter] = newton(U_old, Ep_old, F);     
      
      % test of convergency: 
      if iter==it_max                       % the solver was not succesfull
          d_zeta=max(d_zeta/2,d_zeta_min) ; % decrease of load increment
      else                                  % the solver was succefull
          % update of variables
          step=step+1;
          zeta_old=zeta;
          zeta_hist(step)=zeta;
          settle=-U(2,nv);
          settle_hist(step)=settle;
          if (settle-settle_old)>d_settle      % decrease of the increment
              warning('Too large increment of the settlement.')
              d_zeta=max(d_zeta/2,d_zeta_min) ; 
          end                                   
          settle_old=settle;
          U_old=U;
          Ep_old=Ep;
          disp([' step=', num2str(step), ', settlement=', ...
                num2str(settle),', zeta=', num2str(zeta), ', d_zeta=', ...
                num2str(d_zeta),', iter=', num2str(iter)])
      end
      
      % stopping criteria of the loading process
      if settle_old>=settle_max
          warning('Too large settlement at the point A.')
          break
      end
      
      if d_zeta==d_zeta_min
          warning('Too small load increments.')
          break
      end
      
      if step>=step_max
          warning('Maximal number of steps was achieved.')
          break
      end

  end %true

%
% postprocessing
%

% loading path
figure
hold on;
title('Loading path');
plot(settle_hist(1:(step)), zeta_hist(1:(step)))
hold off

% total displacement
U_total = sqrt(U(1,:).^2 + U(2,:).^2);
figure;
hold on;
title('Total displacement');
patch('Faces',elem','Vertices',coord','FaceVertexCData',U_total',...
      'FaceColor','interp','EdgeColor','none'); 
colorbar;
box on
axis equal
hold off;

% plastic strain 
Ep_norm=sqrt(sum(Ep.*Ep));
Ep_norm_node=transformation(Ep_norm);
figure;
hold on;
title('Plastic strain');
patch('Faces',elem','Vertices',coord','FaceVertexCData',Ep_norm_node',...
      'FaceColor','interp','EdgeColor','none'); 
colorbar;
box on
axis equal
hold off;


  


