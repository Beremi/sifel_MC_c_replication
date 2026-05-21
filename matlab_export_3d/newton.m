% ************************************************************************
%
%  newton  -  an algorithm of the semismooth Newton method
%
% ************************************************************************
%

function [U,        ... % 3*nv displacement vector
          Ep,       ... % 6*nt plastic strains on elements
          iter1,    ... % number of iteration
          last_response ...
          ] = newton(U_init, ...  % 3*nv initial displacement vector
                     Ep_init, ... % 6*nt initial plastic strain
                     F_ext ...    % 3*nv load vector
                    )

  global  volume EU Dir nt nv iMALI jMALI K_elast
  global  Inv_ELAST
  global  tolerance it_max settle_max     % parameters of the Newton method
%
%==========================================================================

% auxilliary notation
  volume6 = ones(6,1)*volume;
  volume36 = ones(36,1)*volume;

% initialization
  U = U_init;        % displacement vector
  E = zeros(6,nt);
  E(:) = EU*U(:);   % strain on elements
  [S,eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,EIG_1,EIG_2,EIG_3,...
    sigma_1,sigma_2,sigma_3,return_type]=constitutive_problem(E,Ep_init);
  Sderiv = stiffness_matrix(eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,...
                            EIG_1,EIG_2,EIG_3,sigma_1,sigma_2,sigma_3);
  F_int = zeros(3,nv);  % vector of internal forces
  dU = zeros(3,nv);     % Newton's increment (in displacement)

% Newton iterations
  it=0;              % iteration number
  while true

      % test on number of iteration
      it=it+1;
      if  it >= it_max
          warning('Maximal number of Newton iteration was achieved.')
          break
      end

      % consistent tangent stiffness matrix
      hMALI = volume36.*Sderiv;
      MALI = sparse( iMALI(:),jMALI(:),hMALI(:) , 6*nt,6*nt );
      K_tangent = EU'*MALI*EU;
      K_tangent = (K_tangent'+K_tangent)/2;

      % vector of internal forces
      F_int(:) = EU'*reshape( volume6.*S , [],1 );

      % Newton's increment
      dU(:) = 0;
      dU(Dir) = K_tangent(Dir,Dir)\(F_ext(Dir)-F_int(Dir));

      % next iteration
      U_new = U + dU;

      % stopping criterion
      q1 = sqrt( abs(dU(:)'*K_elast*dU(:)) );
      q2 = sqrt( abs(  U(:)'*K_elast*U(:)  ) );
      q3 = sqrt( abs( U_new(:)'*K_elast*U_new(:) ) );
      criterion = q1/(q2+q3+eps);

      % update of unknown arrays
      U=U_new;
      E(:) = EU*U(:);
      [S,eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,EIG_1,EIG_2,EIG_3,...
        sigma_1,sigma_2,sigma_3,return_type]=constitutive_problem(E,Ep_init);
      Sderiv = stiffness_matrix(eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,...
                                EIG_1,EIG_2,EIG_3,sigma_1,sigma_2,sigma_3);

      % test on the stopping criterion
      if  criterion < tolerance
          break
      end

      if max(abs(U(:)))>2*settle_max
           warning('Too large displacement in the Newton solver.')
           it=it_max;
           break
      end

  end%  true

  iter1=it;

  Ep = -Inv_ELAST*S;
  Ep(1:6,:)=Ep(1:6,:)+E;

  last_response = struct();
  last_response.E = E;
  last_response.S = S;
  last_response.DS = Sderiv;
  last_response.E_sifel = upstream_to_sifel_rows(E);
  last_response.S_sifel = upstream_to_sifel_rows(S);
  last_response.DS_sifel = upstream_to_sifel_tangent(Sderiv);
  last_response.return_type = return_type;
  last_response.eig_values = [eig_1; eig_2; eig_3];
  last_response.sigma_values = [sigma_1; sigma_2; sigma_3];

end

function A_sifel = upstream_to_sifel_rows(A_upstream)
  upstream_to_sifel = [1, 2, 3, 5, 6, 4];
  A_sifel = A_upstream(upstream_to_sifel, :);
end

function DS_sifel = upstream_to_sifel_tangent(DS_upstream)
  upstream_to_sifel = [1, 2, 3, 5, 6, 4];
  D_upstream = reshape(DS_upstream, 6, 6, []);
  D_sifel = D_upstream(upstream_to_sifel, upstream_to_sifel, :);
  DS_sifel = reshape(D_sifel, 36, []);
end
