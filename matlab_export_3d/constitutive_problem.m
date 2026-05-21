% ************************************************************************
%
%  constitutive_problem   -  solution of the 3D constitutive
%                            Mohr-Coulomb problem
%
% ************************************************************************
%

function [S,...                       % stress tensor 6*nt
          eig_1,eig_2,eig_3,...       % ordered eigenvalues of the trial strain 1*nt
          Eig_1,Eig_2,Eig_3,...       % 1st derivatives ot the eigenvalues 6*nt
          EIG_1,EIG_2,EIG_3,...       % 2nd derivatives ot the eigenvalues 36*nt
          sigma_1,sigma_2,sigma_3,... % ordered principle stresses 1*nt
          return_type]...             % return type used by the export tests
            =        ... %
  constitutive_problem ( ...
            E_new,       ... %  6*nt current strain tensor
            Ep_prev      ... %  6*nt plastic strain from the previous step
                            )
 %
  global nt                                  % number of tetrahedra
  global c_bar sin_phi lame shear            % material parameters
  global iota ELAST IDENT                    % constitutive operators

% .........................................................................

  nt = size(E_new,2);
  sin_psi = sin_phi;                         % associative flow rule

% trial strain
  E_trial = -Ep_prev ;
  E_trial(1:6,:)=E_trial(1:6,:)+E_new;       % trial strain
  E_tr = IDENT*E_trial;                      % trial strain in stress notation
  E_square = ...                             % square of the trial strain
    [ E_tr(1,:).^2         + E_tr(4,:).^2         + E_tr(6,:).^2
      E_tr(2,:).^2         + E_tr(4,:).^2         + E_tr(5,:).^2
      E_tr(3,:).^2         + E_tr(5,:).^2         + E_tr(6,:).^2
      E_tr(1,:).*E_tr(4,:) + E_tr(2,:).*E_tr(4,:) + E_tr(5,:).*E_tr(6,:)
      E_tr(4,:).*E_tr(6,:) + E_tr(2,:).*E_tr(5,:) + E_tr(3,:).*E_tr(5,:)
      E_tr(1,:).*E_tr(6,:) + E_tr(4,:).*E_tr(5,:) + E_tr(3,:).*E_tr(6,:) ];

% ordered eigenvalues of the trial strain
  I1 = E_tr(1,:)+E_tr(2,:)+E_tr(3,:);
  I2 = E_tr(1,:).*E_tr(2,:)+E_tr(1,:).*E_tr(3,:)+...
       E_tr(2,:).*E_tr(3,:)-E_tr(4,:).^2-E_tr(5,:).^2-E_tr(6,:).^2;
  I3 = E_tr(1,:).*E_tr(2,:).*E_tr(3,:)-E_tr(3,:).*E_tr(4,:).^2-...
       E_tr(2,:).*E_tr(6,:).^2-E_tr(1,:).*E_tr(5,:).^2+...
       2*E_tr(4,:).*E_tr(5,:).*E_tr(6,:);
  Q = max(0,(I1.^2-3*I2)/9);
  R = (-2*I1.^3+9*I1.*I2-27*I3)/54;
  test1=(Q==0);
  theta0 = zeros(1,nt);
  theta0(~test1)=R(~test1)./sqrt(Q(~test1).^3);
  theta=acos(min(max(theta0,-1),1))/3;

  eig_1 = -2*sqrt(Q).*cos(theta+2*pi/3)+I1/3;
  eig_2 = -2*sqrt(Q).*cos(theta-2*pi/3)+I1/3;
  eig_3 = -2*sqrt(Q).*cos(theta)+I1/3;

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
  denom_a=4*(lame+2*shear/3)*sin_phi*sin_psi;

  lambda_s=f_tr/denom_s ;
  lambda_l=(shear*((1+sin_phi)*(eig_1+eig_2)-2*(1-sin_phi)*eig_3)+ ...
            2*lame*sin_phi*trace_E-c_bar)/denom_l ;
  lambda_r=(shear*(2*(1+sin_phi)*eig_1-(1-sin_phi)*(eig_2+eig_3))+ ...
            2*lame*sin_phi*trace_E-c_bar)/denom_r ;
  lambda_a=(2*(lame+2*shear/3)*sin_phi*trace_E-c_bar)/denom_a ;

% initialization of unknowns variables
  S = zeros(6,nt);
  sigma_1 = zeros(1,nt);
  sigma_2 = zeros(1,nt);
  sigma_3 = zeros(1,nt);   % principal stresses
  Eig_1 = zeros(6,nt);
  Eig_2 = zeros(6,nt);
  Eig_3 = zeros(6,nt);     % eigenprojections of the trial strain (1st der.)
  EIG_1 = zeros(36,nt);
  EIG_2 = zeros(36,nt);
  EIG_3 = zeros(36,nt);    % the second derivatives of the eigenvalues
  return_type = zeros(1,nt);

% auxilliary array for the second derivatives of the eigenprojections
  DER_E_square = ...       % 36*nt array (fourth order tensors)
    [ 2*E_tr(1,:)   ; zeros(1,nt); zeros(1,nt); E_tr(4,:)                ; zeros(1,nt)           ; E_tr(6,:)
      zeros(1,nt); 2*E_tr(2,:)   ; zeros(1,nt); E_tr(4,:)                ; E_tr(5,:)            ; zeros(1,nt)
      zeros(1,nt); zeros(1,nt); 2*E_tr(3,:)   ; zeros(1,nt)              ; E_tr(5,:)            ; E_tr(6,:)
      E_tr(4,:)  ; E_tr(4,:)     ; zeros(1,nt); 0.5*(E_tr(1,:)+E_tr(2,:)); 0.5*E_tr(6,:)        ; 0.5*E_tr(5,:)
      zeros(1,nt); E_tr(5,:)     ; E_tr(5,:) ; 0.5*E_tr(6,:)             ; 0.5*(E_tr(2,:)+E_tr(3,:)); 0.5*E_tr(4,:)
      E_tr(6,:)  ; zeros(1,nt)   ; E_tr(6,:) ; 0.5*E_tr(5,:)             ; 0.5*E_tr(4,:)        ; 0.5*(E_tr(1,:)+E_tr(3,:)) ];

% elastic response
  test_el=(f_tr<=0);
  nt_el=length(sigma_1(test_el));
  sigma_1(test_el)=lame*trace_E(test_el)+2*shear*eig_1(test_el);
  sigma_2(test_el)=lame*trace_E(test_el)+2*shear*eig_2(test_el);
  sigma_3(test_el)=lame*trace_E(test_el)+2*shear*eig_3(test_el);
  S(:,test_el) = ELAST*E_trial(:,test_el) ;

% return to the smooth portion of the yield surface
  test_s=(lambda_s<=min(gamma_sl,gamma_sr))&(~test_el);
  lambda_s=lambda_s(test_s);
  nt_s=length(sigma_1(test_s));
  eig_1_s=eig_1(test_s);
  eig_2_s=eig_2(test_s);
  eig_3_s=eig_3(test_s);
  E_square_s=E_square(:,test_s);
  E_tr_s=E_tr(:,test_s);
  denom_s1=(eig_1_s-eig_2_s).*(eig_1_s-eig_3_s);
  denom_s2=(eig_2_s-eig_1_s).*(eig_2_s-eig_3_s);
  denom_s3=(eig_3_s-eig_1_s).*(eig_3_s-eig_2_s);
  Eig_1_s=(ones(6,1)*(1./denom_s1)).*(E_square_s-...
          (ones(6,1)*(eig_2_s+eig_3_s)).*E_tr_s+iota*(eig_2_s.*eig_3_s));
  Eig_2_s=(ones(6,1)*(1./denom_s2)).*(E_square_s-...
          (ones(6,1)*(eig_1_s+eig_3_s)).*E_tr_s+iota*(eig_1_s.*eig_3_s));
  Eig_3_s=(ones(6,1)*(1./denom_s3)).*(E_square_s-...
          (ones(6,1)*(eig_1_s+eig_2_s)).*E_tr_s+iota*(eig_1_s.*eig_2_s));
  Eig_1(:,test_s)=Eig_1_s;
  Eig_2(:,test_s)=Eig_2_s;
  Eig_3(:,test_s)=Eig_3_s;
  [EIG_1(:,test_s),EIG_2(:,test_s),EIG_3(:,test_s)] = ...
      eigenprojection_derivatives_smooth(DER_E_square(:,test_s),IDENT,...
        Eig_1_s,Eig_2_s,Eig_3_s,eig_1_s,eig_2_s,eig_3_s,...
        denom_s1,denom_s2,denom_s3);
  sigma_1(test_s)=lame*trace_E(test_s)+2*shear*eig_1(test_s)-...
                  lambda_s*(2*lame*sin_psi+2*shear*(1+sin_psi));
  sigma_2(test_s)=lame*trace_E(test_s)+2*shear*eig_2(test_s)-...
                  lambda_s*(2*lame*sin_psi);
  sigma_3(test_s)=lame*trace_E(test_s)+2*shear*eig_3(test_s)-...
                  lambda_s*(2*lame*sin_psi-2*shear*(1-sin_psi));
  S(:,test_s)=(ones(6,1)*sigma_1(test_s)).*Eig_1(:,test_s)+...
              (ones(6,1)*sigma_2(test_s)).*Eig_2(:,test_s)+...
              (ones(6,1)*sigma_3(test_s)).*Eig_3(:,test_s);
  return_type(test_s)=1;

% return to the left edge of the yield surface
  test_l=(gamma_sl<gamma_sr)&(lambda_l>=gamma_sl)&(lambda_l<=gamma_la)&...
         (~(test_el|test_s));
  lambda_l=lambda_l(test_l);
  nt_l=length(sigma_1(test_l));
  eig_1_l=eig_1(test_l);
  eig_2_l=eig_2(test_l);
  eig_3_l=eig_3(test_l);
  E_square_l=E_square(:,test_l);
  E_tr_l=E_tr(:,test_l);
  denom_l3=(eig_3_l-eig_1_l).*(eig_3_l-eig_2_l);
  Eig_3_l=(ones(6,1)*(1./denom_l3)).*(E_square_l-...
          (ones(6,1)*(eig_1_l+eig_2_l)).*E_tr_l+iota*(eig_1_l.*eig_2_l));
  Eig_12_l=[ones(3,nt_l); zeros(3,nt_l)]-Eig_3_l;
  Eig_3(:,test_l)=Eig_3_l;
  EIG_3(:,test_l)=eigenprojection_derivative_left(DER_E_square(:,test_l),...
    IDENT,E_tr_l,Eig_12_l,Eig_3_l,eig_1_l,eig_2_l,eig_3_l,denom_l3);
  sigma_1(test_l)=lame*trace_E(test_l)+...
                  shear*(eig_1(test_l)+eig_2(test_l))-...
                  lambda_l*(2*lame*sin_psi+shear*(1+sin_psi));
  sigma_2(test_l)=sigma_1(test_l);
  sigma_3(test_l)=lame*trace_E(test_l)+2*shear*eig_3(test_l)-...
                  lambda_l*(2*lame*sin_psi-2*shear*(1-sin_psi));
  S(:,test_l)=(ones(6,1)*sigma_1(test_l)).*Eig_12_l+...
              (ones(6,1)*sigma_3(test_l)).*Eig_3(:,test_l);
  return_type(test_l)=2;

% return to the right edge of the yield surface
  test_r=(gamma_sl>gamma_sr)&(lambda_r>=gamma_sr)&(lambda_r<=gamma_ra)&...
         (~(test_el|test_s));
  lambda_r=lambda_r(test_r);
  nt_r=length(sigma_1(test_r));
  eig_1_r=eig_1(test_r);
  eig_2_r=eig_2(test_r);
  eig_3_r=eig_3(test_r);
  E_square_r=E_square(:,test_r);
  E_tr_r=E_tr(:,test_r);
  denom_r1=(eig_1_r-eig_2_r).*(eig_1_r-eig_3_r);
  Eig_1_r=(ones(6,1)*(1./denom_r1)).*(E_square_r-...
          (ones(6,1)*(eig_2_r+eig_3_r)).*E_tr_r+iota*(eig_2_r.*eig_3_r));
  Eig_23_r=[ones(3,nt_r); zeros(3,nt_r)]-Eig_1_r;
  Eig_1(:,test_r)=Eig_1_r;
  EIG_1(:,test_r)=eigenprojection_derivative_right(DER_E_square(:,test_r),...
    IDENT,E_tr_r,Eig_1_r,Eig_23_r,eig_1_r,eig_2_r,eig_3_r,denom_r1);
  sigma_1(test_r)=lame*trace_E(test_r)+2*shear*eig_1(test_r)-...
                  lambda_r*(2*lame*sin_psi+2*shear*(1+sin_psi));
  sigma_3(test_r)=lame*trace_E(test_r)+...
                  shear*(eig_2(test_r)+eig_3(test_r))-...
                  lambda_r*(2*lame*sin_psi-shear*(1-sin_psi));
  sigma_2(test_r)=sigma_3(test_r);
  S(:,test_r)=(ones(6,1)*sigma_1(test_r)).*Eig_1(:,test_r)+...
              (ones(6,1)*sigma_3(test_r)).*Eig_23_r;
  return_type(test_r)=3;

% return to the apex of the yield surface
  test_a=~(test_el|test_s|test_l|test_r);
  lambda_a=lambda_a(test_a);
  nt_a=length(lambda_a);
  if (nt~=nt_el+nt_s+nt_l+nt_r+nt_a), warning('number of elements!'), end
  sigma_1(test_a)=ones(1,nt_a)*c_bar/(2*sin_phi);
  sigma_2(test_a)=sigma_1(test_a);
  sigma_3(test_a)=sigma_1(test_a);
  S(:,test_a)=iota*sigma_1(test_a);
  return_type(test_a)=4;

end

function [EIG_1,EIG_2,EIG_3] = eigenprojection_derivatives_smooth( ...
  DER_E_square,IDENT,Eig_1,Eig_2,Eig_3,eig_1,eig_2,eig_3,...
  denom_1,denom_2,denom_3)

  E1_x_E1 = column_outer(Eig_1,Eig_1);
  E2_x_E2 = column_outer(Eig_2,Eig_2);
  E3_x_E3 = column_outer(Eig_3,Eig_3);

  EIG_1 = (ones(36,1)*(1./denom_1)).*(DER_E_square-...
          IDENT(:)*(eig_2+eig_3)-...
          (ones(36,1)*(2*eig_1-eig_2-eig_3)).*E1_x_E1-...
          (ones(36,1)*(eig_2-eig_3)).*(E2_x_E2-E3_x_E3));
  EIG_2 = (ones(36,1)*(1./denom_2)).*(DER_E_square-...
          IDENT(:)*(eig_1+eig_3)-...
          (ones(36,1)*(2*eig_2-eig_1-eig_3)).*E2_x_E2-...
          (ones(36,1)*(eig_1-eig_3)).*(E1_x_E1-E3_x_E3));
  EIG_3 = (ones(36,1)*(1./denom_3)).*(DER_E_square-...
          IDENT(:)*(eig_1+eig_2)-...
          (ones(36,1)*(2*eig_3-eig_1-eig_2)).*E3_x_E3-...
          (ones(36,1)*(eig_1-eig_2)).*(E1_x_E1-E2_x_E2));

end

function EIG_3 = eigenprojection_derivative_left( ...
  DER_E_square,IDENT,E_tr,Eig_12,Eig_3,eig_1,eig_2,eig_3,denom_3)

  E3_x_E3 = column_outer(Eig_3,Eig_3);
  E12_x_E12 = column_outer(Eig_12,Eig_12);
  E12_x_Etr = column_outer(Eig_12,E_tr);
  Etr_x_E12 = column_outer(E_tr,Eig_12);
  E12_x_E3 = column_outer(Eig_12,Eig_3);
  E3_x_E12 = column_outer(Eig_3,Eig_12);

  EIG_3 = (ones(36,1)*(1./denom_3)).*(DER_E_square-...
          IDENT(:)*(eig_1+eig_2)-(Etr_x_E12+E12_x_Etr)+...
          (ones(36,1)*(eig_1+eig_2)).*E12_x_E12+...
          (ones(36,1)*(eig_1+eig_2-2*eig_3)).*E3_x_E3+...
          (ones(36,1)*eig_3).*(E12_x_E3+E3_x_E12));

end

function EIG_1 = eigenprojection_derivative_right( ...
  DER_E_square,IDENT,E_tr,Eig_1,Eig_23,eig_1,eig_2,eig_3,denom_1)

  E1_x_E1 = column_outer(Eig_1,Eig_1);
  E23_x_E23 = column_outer(Eig_23,Eig_23);
  E23_x_Etr = column_outer(Eig_23,E_tr);
  Etr_x_E23 = column_outer(E_tr,Eig_23);
  E23_x_E1 = column_outer(Eig_23,Eig_1);
  E1_x_E23 = column_outer(Eig_1,Eig_23);

  EIG_1 = (ones(36,1)*(1./denom_1)).*(DER_E_square-...
          IDENT(:)*(eig_2+eig_3)-(Etr_x_E23+E23_x_Etr)+...
          (ones(36,1)*(eig_2+eig_3)).*E23_x_E23+...
          (ones(36,1)*(eig_2+eig_3-2*eig_1)).*E1_x_E1+...
          (ones(36,1)*eig_1).*(E23_x_E1+E1_x_E23));

end

function M = column_outer(A,B)

  M = reshape(reshape(A,6,1,[]).*reshape(B,1,6,[]),36,[]);

end
