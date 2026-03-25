% ************************************************************************
%
%  constitutive_problem   -  solution of the constitutive 
%                            Mohr-Coulomb problem 
%
% ************************************************************************
%                          

function [S,...                       % stress tensor 4xnt 
          eig_1,eig_2,eig_3,...       % ordered eigenvalues of the trial strain 1xnt
          Eig_1,Eig_2,Eig_3,...       % 1st derivatives ot the eigenvalues 4xnt
          EIG_1,EIG_2,EIG_3,...       % 2nd derivatives ot the eigenvalues 9xnt
          sigma_1,sigma_2,sigma_3]... % ordered principle stresses 1xnt
            =        ... %           
  constitutive_problem ( ...
            E_new,       ... %  3*nt current strain tensor
            Ep_prev      ... %  4*nt plastic strain from the previous step  
                            )
 %
  global nt                                  % number of triangles
  global c_bar sin_phi lame shear            % material parameters  
  global iota ELAST                          % constitutive operators

% .........................................................................

% trial strain
  E_trial = -Ep_prev ;
  E_trial(1:3,:)=E_trial(1:3,:)+E_new;      % trial strain
  
% inordered eigenvalues of the trial strain
  I1=E_trial(1,:)+E_trial(2,:);
  I2=sqrt((E_trial(1,:)-E_trial(2,:)).^2+E_trial(3,:).^2);
  eig0_1=(I1+I2)/2 ;
  eig0_2=(I1-I2)/2 ;
  eig0_3=E_trial(4,:) ;
  test1=(I2==0); % first two eigenvalues are the same

% inordered first derivatives of the trial stress
  Eig0_1 = zeros(4,nt);
  Eig0_1(1:3,~test1) = [E_trial(1,~test1)-eig0_2(~test1)
                        E_trial(2,~test1)-eig0_2(~test1)
                        E_trial(3,~test1)/2]./...
                       (ones(3,1)*I2(~test1));
  Eig0_1(1:2,test1) = ones(2,length(eig0_1(test1)));               
  Eig0_2 = [ones(2,nt); zeros(2,nt)]-Eig0_1;
  Eig0_3 = [zeros(3,nt) ; ones(1,nt)];

% inordered second derivatives of the trial stress
  EIG0_1 = zeros(9,nt);
  EIG0_2 = zeros(9,nt);
  EIG0_3 = zeros(9,nt);   
  EIG0_1(:,~test1)=...
  [1-Eig0_1(1,~test1).*Eig0_1(1,~test1)-Eig0_2(1,~test1).*Eig0_2(1,~test1)
    -Eig0_1(2,~test1).*Eig0_1(1,~test1)-Eig0_2(2,~test1).*Eig0_2(1,~test1)
    -Eig0_1(3,~test1).*Eig0_1(1,~test1)-Eig0_2(3,~test1).*Eig0_2(1,~test1)
    -Eig0_1(1,~test1).*Eig0_1(2,~test1)-Eig0_2(1,~test1).*Eig0_2(2,~test1)
   1-Eig0_1(2,~test1).*Eig0_1(2,~test1)-Eig0_2(2,~test1).*Eig0_2(2,~test1)
    -Eig0_1(3,~test1).*Eig0_1(2,~test1)-Eig0_2(3,~test1).*Eig0_2(2,~test1)
    -Eig0_1(1,~test1).*Eig0_1(3,~test1)-Eig0_2(1,~test1).*Eig0_2(3,~test1)
    -Eig0_1(2,~test1).*Eig0_1(3,~test1)-Eig0_2(2,~test1).*Eig0_2(3,~test1)
  1/2-Eig0_1(3,~test1).*Eig0_1(3,~test1)-Eig0_2(3,~test1).*Eig0_2(3,~test1)...
  ]./(ones(9,1)*I2(~test1));
  EIG0_2(:,~test1)=-EIG0_1(:,~test1);
 
% reordering of eigenvalues and their derivatives  
  
  eig_1 = eig0_1; 
  eig_2 = eig0_2;   
  eig_3 = eig0_3;   % ordered eigenvalues of the trial strain
  Eig_1 = Eig0_1;
  Eig_2 = Eig0_2;
  Eig_3 = Eig0_3;   % eigenprojections of the trial strain (1st der.)
  EIG_1 = EIG0_1;
  EIG_2 = EIG0_2;
  EIG_3 = EIG0_3;   % the second derivatives of the eigenvalues
       
  test2=(eig0_1>=eig0_3)&(eig0_3>eig0_2);
  eig_2(test2)=eig0_3(test2);
  eig_3(test2)=eig0_2(test2);
  Eig_2(:,test2)=Eig0_3(:,test2);
  Eig_3(:,test2)=Eig0_2(:,test2);
  EIG_2(:,test2)=EIG0_3(:,test2);
  EIG_3(:,test2)=EIG0_2(:,test2);
   
  test3=(eig0_3>eig0_1);
  eig_1(test3)=eig0_3(test3);
  eig_2(test3)=eig0_1(test3);
  eig_3(test3)=eig0_2(test3);   % reordered eigenvalues of the trial strain
  Eig_1(:,test3)=Eig0_3(:,test3);
  Eig_2(:,test3)=Eig0_1(:,test3);
  Eig_3(:,test3)=Eig0_2(:,test3);
  EIG_1(:,test3)=EIG0_3(:,test3);
  EIG_2(:,test3)=EIG0_1(:,test3);
  EIG_3(:,test3)=EIG0_2(:,test3); % reordered derivatives
  
% critical values 
  trace_E=eig_1+eig_2+eig_3;        % trace of E_trial
  f_tr=2*shear*((1+sin_phi)*eig_1-(1-sin_phi)*eig_3)+ ...
       2*lame*sin_phi*trace_E-c_bar; % test on admissibility
  gamma_sl=(eig_1-eig_2)/(1+sin_phi);
  gamma_sr=(eig_2-eig_3)/(1-sin_phi);
  gamma_la=(eig_1+eig_2-2*eig_3)/(3-sin_phi);
  gamma_ra=(2*eig_1-eig_2-eig_3)/(3+sin_phi);
  
% candidates on plastic multipliers
  denom_s=4*lame*sin_phi^2+2*shear*(1+sin_phi)^2+2*shear*(1-sin_phi)^2;
  denom_l=4*lame*sin_phi^2+  shear*(1+sin_phi)^2+2*shear*(1-sin_phi)^2;
  denom_r=4*lame*sin_phi^2+2*shear*(1+sin_phi)^2+  shear*(1-sin_phi)^2;
  
  lambda_s=f_tr/denom_s ;
  lambda_l=(shear*((1+sin_phi)*(eig_1+eig_2)-2*(1-sin_phi)*eig_3)+ ...
            2*lame*sin_phi*trace_E-c_bar)/denom_l ;
  lambda_r=(shear*(2*(1+sin_phi)*eig_1-(1-sin_phi)*(eig_2+eig_3))+ ...
            2*lame*sin_phi*trace_E-c_bar)/denom_r ;
  
% initialization of unknowns variables
  S = zeros(4,nt);
  sigma_1 = zeros(1,nt); 
  sigma_2 = zeros(1,nt);   
  sigma_3 = zeros(1,nt);   % principal stresses
  
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
  sigma_1(test_s)=lame*trace_E(test_s)+2*shear*eig_1(test_s)-...
                  lambda_s*(2*lame*sin_phi+2*shear*(1+sin_phi));
  sigma_2(test_s)=lame*trace_E(test_s)+2*shear*eig_2(test_s)-...
                  lambda_s*(2*lame*sin_phi);
  sigma_3(test_s)=lame*trace_E(test_s)+2*shear*eig_3(test_s)-...
                  lambda_s*(2*lame*sin_phi-2*shear*(1-sin_phi));
  S(:,test_s)=(ones(4,1)*sigma_1(test_s)).*Eig_1(:,test_s)+...
              (ones(4,1)*sigma_2(test_s)).*Eig_2(:,test_s)+...
              (ones(4,1)*sigma_3(test_s)).*Eig_3(:,test_s);
              
% return to the left edge of the yield surface             
  test_l=(gamma_sl<gamma_sr)&(lambda_l>=gamma_sl)&(lambda_l<=gamma_la)&...
         (~(test_el|test_s));
  lambda_l=lambda_l(test_l);  
  nt_l=length(sigma_1(test_l));
  sigma_1(test_l)=lame*trace_E(test_l)+...
                  shear*(eig_1(test_l)+eig_2(test_l))-...
                  lambda_l*(2*lame*sin_phi+shear*(1+sin_phi));
  sigma_2(test_l)=sigma_1(test_l);
  sigma_3(test_l)=lame*trace_E(test_l)+2*shear*eig_3(test_l)-...
                  lambda_l*(2*lame*sin_phi-2*shear*(1-sin_phi));
  S(:,test_l)=(ones(4,1)*sigma_1(test_l)).*...
                                (Eig_1(:,test_l)+Eig_2(:,test_l))+...
              (ones(4,1)*sigma_3(test_l)).*Eig_3(:,test_l);   
          
% return to the right edge of the yield surface             
  test_r=(gamma_sl>gamma_sr)&(lambda_r>=gamma_sr)&(lambda_r<=gamma_ra)&...
         (~(test_el|test_s));
  lambda_r=lambda_r(test_r);     
  nt_r=length(sigma_1(test_r));
  sigma_1(test_r)=lame*trace_E(test_r)+2*shear*eig_1(test_r)-...
                  lambda_r*(2*lame*sin_phi+2*shear*(1+sin_phi));
  sigma_3(test_r)=lame*trace_E(test_r)+...
                  shear*(eig_2(test_r)+eig_3(test_r))-...
                  lambda_r*(2*lame*sin_phi-shear*(1-sin_phi));
  sigma_2(test_r)=sigma_3(test_r);
  S(:,test_r)=(ones(4,1)*sigma_1(test_r)).*Eig_1(:,test_r)+...                                
              (ones(4,1)*sigma_3(test_r)).*...
                  (Eig_2(:,test_r)+Eig_3(:,test_r));        
              
% return to the apex of the yield surface
  test_a=~(test_el|test_s|test_l|test_r);
  nt_a=length(sigma_1(test_a));
  if (nt~=nt_el+nt_s+nt_l+nt_r+nt_a), warning('number of elements!'), end
  sigma_1(test_a)=ones(1,nt_a)*c_bar/(2*sin_phi);
  sigma_2(test_a)=sigma_1(test_a);
  sigma_3(test_a)=sigma_1(test_a);      
  S(:,test_a)=iota*sigma_1(test_a);  
  
 end
