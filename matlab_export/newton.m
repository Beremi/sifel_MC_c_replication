% ************************************************************************
%
%  newton  -  an algorithm of the semismooth Newton method
%
% ************************************************************************
%

function [U,        ... % 2*nv displacement vector
          Ep,       ... % 4*nt plastic strains on elements
          iter1     ... % number of iteration
          ] = newton(U_init, ...  % 2*nv initial displacement vector
                     Ep_init, ... % 4*nt initial plastic strain
                     F_ext ...    % 2*nv load vector
                    )                   
      
  global  shear bulk lame c_bar sin_phi           % material parameters
  global  area EU Dir nt nv iMALI jMALI K_elast
  global  DEV Dev VOL iota ELAST Elast Inv_ELAST  % constitutive operators 
  global  IDENT Ident
  global  tolerance it_max settle_max     % parameters of the Newton method
%
%==========================================================================
  
% auxilliary notation
  area3 = ones(3,1)*area ;
  area9 = ones(9,1)*area ;

% initialization
  U = U_init ;       % displacement vector
  E = zeros(3,nt) ;
  E(:) = EU*U(:) ;   % strain on elements
  [S,eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,EIG_1,EIG_2,EIG_3,...
    sigma_1,sigma_2,sigma_3]=constitutive_problem (E,Ep_init) ; 
  Sderiv = stiffness_matrix(eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,...
                            EIG_1,EIG_2,EIG_3,sigma_1,sigma_2,sigma_3);
  F_int = zeros(2,nv) ;  % vector of internal forces
  dU = zeros(2,nv) ; % Newton's increment (in displacement)
 
% Newton iterations
  it=0;              % iteration number
  while true   
      
      % test on number of iteration
      it=it+1; 
      if  it >= it_max
          % plastic strain 
          Ep=-Inv_ELAST*S;
          Ep(1:3,:)=Ep(1:3,:)+E;
          warning('Maximal number of Newton iteration was achieved.')
          break
      end       
      
      % consistent tangent stiffness matrix
      hMALI = area9.*Sderiv ; 
      MALI = sparse( iMALI(:),jMALI(:),hMALI(:) , 3*nt,3*nt ) ;
      K_tangent = EU'*MALI*EU ;   
      K_tangent = (K_tangent'+K_tangent)/2  ;   
 
      % vector of internal forces
      F_int(:) = EU'*reshape( area3.*S(1:3,:) , [],1 ) ;
      %norm(F_ext(Dir)-F_int(Dir))/norm(F_ext(Dir))
      
      % Newton's increment
      dU(Dir) = K_tangent(Dir,Dir)\(F_ext(Dir)-F_int(Dir)); 
            
      % next iteration
      U_new = U + dU ;

      % stopping criterion 
      q1 = sqrt( dU(:)'*K_elast*dU(:) ) ;
      q2 = sqrt(  U(:)'*K_elast*U(:)  ) ;
      q3 = sqrt( U_new(:)'*K_elast*U_new(:) ) ;
      criterion = q1/(q2+q3);
      
      % update of unknown arrays
      U=U_new;
      E(:) = EU*U(:) ;
      [S,eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,EIG_1,EIG_2,EIG_3,...
        sigma_1,sigma_2,sigma_3]=constitutive_problem (E,Ep_init) ; 
      Sderiv = stiffness_matrix(eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,...
                            EIG_1,EIG_2,EIG_3,sigma_1,sigma_2,sigma_3);
      
      % test on the stopping criterion
      if  criterion < tolerance
          % plastic strain 
          Ep=-Inv_ELAST*S;
          Ep(1:3,:)=Ep(1:3,:)+E;
          break
      end      
   
      if abs(U(2,nv))>2*settle_max
           % plastic strain 
           Ep=-Inv_ELAST*S;
           Ep(1:3,:)=Ep(1:3,:)+E;
           warning('Too large displacement in the Newton solver.')
           it=it_max;
           break
      end

  end%  true
  
  iter1=it;  
end 