% ************************************************************************
%
%  stiffness_matrix   -  assembling of the tangent stiffness matrix for the
%                              3D Mohr-Coulomb problem
%
% ************************************************************************

function   Sderiv          ... %  36*nt consistent tangent stiffness matrix
             =             ... %
  stiffness_matrix...
     (eig_1,eig_2,eig_3,...    % ordered eigenvalues of the trial strain 1*nt
      Eig_1,Eig_2,Eig_3,...    % 1st derivatives ot the eigenvalues 6*nt
      EIG_1,EIG_2,EIG_3,....   % 2nd derivatives ot the eigenvalues 36*nt
      sigma_1,sigma_2,sigma_3) % ordered principle stresses 1*nt
 %
  global nt                               % number of tetrahedra
  global c_bar sin_phi lame shear         % material parameters
  global iota Elast                       % elastic stiffness matrix

% .........................................................................

  sin_psi = sin_phi;                      % associative flow rule

% critical values
  trace_E=eig_1+eig_2+eig_3;        % trace of E_trial
  f_tr=2*shear*((1+sin_phi)*eig_1-(1-sin_phi)*eig_3)+ ...
       2*lame*sin_phi*trace_E-c_bar; % test on admissibility
  gamma_sl=(eig_1-eig_2)/(1+sin_psi);
  gamma_sr=(eig_2-eig_3)/(1-sin_psi);
  gamma_la=(eig_1+eig_2-2*eig_3)/(3-sin_psi);
  gamma_ra=(2*eig_1-eig_2-eig_3)/(3+sin_psi);

% candidates on plastic multipliers
  denom_s=4*lame*sin_phi*sin_psi+4*shear*(1+sin_phi*sin_psi);
  denom_l=4*lame*sin_phi*sin_psi+shear*(1+sin_phi)*(1+sin_psi)+...
          2*shear*(1-sin_phi)*(1-sin_psi);
  denom_r=4*lame*sin_phi*sin_psi+2*shear*(1+sin_phi)*(1+sin_psi)+...
          shear*(1-sin_phi)*(1-sin_psi);

  lambda_s=f_tr/denom_s ;
  lambda_l=(shear*((1+sin_phi)*(eig_1+eig_2)-2*(1-sin_phi)*eig_3)+ ...
            2*lame*sin_phi*trace_E-c_bar)/denom_l ;
  lambda_r=(shear*(2*(1+sin_phi)*eig_1-(1-sin_phi)*(eig_2+eig_3))+ ...
            2*lame*sin_phi*trace_E-c_bar)/denom_r ;

% initialization
  Sderiv = zeros(36,nt);

% elastic response
  test_el=(f_tr<=0);
  nt_el=length(sigma_1(test_el));
  Sderiv(:,test_el)=Elast(:)*ones(1,nt_el);

% return to the smooth portion of the yield surface
  test_s=(lambda_s<=min(gamma_sl,gamma_sr))&(~test_el);
  nt_s=length(sigma_1(test_s));
  mat1_s=(ones(36,1)*sigma_1(test_s)).*EIG_1(:,test_s)+...
         (ones(36,1)*sigma_2(test_s)).*EIG_2(:,test_s)+...
         (ones(36,1)*sigma_3(test_s)).*EIG_3(:,test_s);
  mat2_s=lame*column_outer(iota*ones(1,nt_s),iota*ones(1,nt_s));
  mat3_s=2*shear*column_outer(Eig_1(:,test_s),Eig_1(:,test_s));
  mat4_s=2*shear*column_outer(Eig_2(:,test_s),Eig_2(:,test_s));
  mat5_s=2*shear*column_outer(Eig_3(:,test_s),Eig_3(:,test_s));
  Eig_6=2*shear*((1+sin_phi)*Eig_1(:,test_s)-...
                 (1-sin_phi)*Eig_3(:,test_s))+...
        2*lame*sin_phi*iota*ones(1,nt_s);
  mat6_s=column_outer(Eig_6,Eig_6)/denom_s;
  Sderiv(:,test_s)=mat1_s+mat2_s+mat3_s+mat4_s+mat5_s-mat6_s;

% return to the left edge of the yield surface
  test_l=(gamma_sl<gamma_sr)&(lambda_l>=gamma_sl)&(lambda_l<=gamma_la)&...
         (~(test_el|test_s));
  nt_l=length(sigma_1(test_l));
  Eig_12=[ones(3,nt_l); zeros(3,nt_l)]-Eig_3(:,test_l);
  mat1_l=(ones(36,1)*(sigma_3(test_l)-sigma_1(test_l))).*EIG_3(:,test_l);
  mat2_l=lame*column_outer(iota*ones(1,nt_l),iota*ones(1,nt_l));
  mat3_l=shear*column_outer(Eig_12,Eig_12);
  mat5_l=2*shear*column_outer(Eig_3(:,test_l),Eig_3(:,test_l));
  Eig_6=shear*((1+sin_phi)*Eig_12-2*(1-sin_phi)*Eig_3(:,test_l))+...
        2*lame*sin_phi*iota*ones(1,nt_l);
  mat6_l=column_outer(Eig_6,Eig_6)/denom_l;
  Sderiv(:,test_l)=mat1_l+mat2_l+mat3_l+mat5_l-mat6_l;

% return to the right edge of the yield surface
  test_r=(gamma_sl>gamma_sr)&(lambda_r>=gamma_sr)&(lambda_r<=gamma_ra)&...
         (~(test_el|test_s));
  nt_r=length(sigma_1(test_r));
  Eig_23=[ones(3,nt_r); zeros(3,nt_r)]-Eig_1(:,test_r);
  mat1_r=(ones(36,1)*(sigma_1(test_r)-sigma_3(test_r))).*EIG_1(:,test_r);
  mat2_r=lame*column_outer(iota*ones(1,nt_r),iota*ones(1,nt_r));
  mat3_r=2*shear*column_outer(Eig_1(:,test_r),Eig_1(:,test_r));
  mat5_r=shear*column_outer(Eig_23,Eig_23);
  Eig_6=shear*(2*(1+sin_phi)*Eig_1(:,test_r)-(1-sin_phi)*Eig_23)+...
        2*lame*sin_phi*iota*ones(1,nt_r);
  mat6_r=column_outer(Eig_6,Eig_6)/denom_r;
  Sderiv(:,test_r)=mat1_r+mat2_r+mat3_r+mat5_r-mat6_r;

% return to the apex of the yield surface
  test_a=~(test_el|test_s|test_l|test_r);
  nt_a=length(sigma_1(test_a));
  if (nt~=nt_el+nt_s+nt_l+nt_r+nt_a), warning('number of elements!'), end
  Sderiv(:,test_a)=zeros(36,nt_a);

end

function M = column_outer(A,B)

  M = reshape(reshape(A,6,1,[]).*reshape(B,1,6,[]),36,[]);

end
