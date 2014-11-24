
% QLM

% Calculating Equilibria, Stability Test, NPZ Model QLM 
% ND = 1.5 or 8 or 16; 2 Parameter sets (Exp.II,III)
% Diploma Thesis P.Fuhs (CAU)
% Matlab Student Version V 7.10.0 (R2010a)

clear all;
clc;
clf;

delete('QLM_diary');
diary on; 
diary('QLM_diary');

countNotifications = 0;

disp 'Calculation for QLM';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = NaN(1,400);         % Layer depth from 1 to 400 m, grid = 1m 

for j= 1:400
    M(j)=j;
end; 

% All units as in thesis.
% All domains as in thesis.

% Vector ND values

ND_vector = [1.5, 8, 16];

% Parameter sets: Exp. II, III

mu_vector = [1.679, 3.775];
Phi_P_vector = [0.033, 0.017];
k_vector = [0.636, 0.560];
g_vector = [0.982, 1.560];
eps_vector = [3.473, 3.740];
Phi_Z_vector = [0.010, 0.011];
beta_vector = [1.003, 1.045];
m_r_vector = [0.107, 1.270];
Omega_vector = [0.990, 1.006];
gamma_vector = [0.990, 0.870];
Phi_P_opt_vector = [0.168, 0.122];
Phi_Z_opt_vector = [0.257, 0.098];


% Create tables with results:
% 
% 6 tables: ND = 1.5,8,16; Exp. II,III => 6 cases
% x 2 for QLM 2 and QLM 3
% table names: ExpNumber_NDValue_QLMNumber_table, e.g.II_8_QLM3_table
% 1 in table names refers to ND = 1.5
% lines in result matrices:
posM = 1;           % line 1: M
posND = 2;          % line 2: ND
posExp = 3;         % line 3: Exp
posN = 4;           % line 4: N
posdotN = 5;        % line 5: dotN
posP = 6;           % line 6: P
posdotP = 7;        % line 7: dotP
posZ = 8;           % line 8: Z
posdotZ = 9;        % line 9: dotZ
posG1 = 10;         % line 10: G_1
posxi1 = 11;        % line 11: xi_1
posxi3 = 12;        % line 12: xi_3
posxi4 = 13;        % line 13: xi_4
posxi5 = 14;        % line 14: xi_5
posxi6 = 15;        % line 15: xi_6
posxi7 = 16;        % line 16: xi_7
posu = 17;          % line 17: u
posG = 18;          % line 18: G
posFN = 19;         % line 19: FN
posStable = 20;     % line 20: Stability

II_1_QLM2_table = NaN(20,400);    II_1_QLM3_table = NaN(20,400); 
II_8_QLM2_table = NaN(20,400);    II_8_QLM3_table = NaN(20,400); 
II_16_QLM2_table = NaN(20,400);   II_16_QLM3_table = NaN(20,400); 
III_1_QLM2_table = NaN(20,400);   III_1_QLM3_table = NaN(20,400); 
III_8_QLM2_table = NaN(20,400);   III_8_QLM3_table = NaN(20,400);  
III_16_QLM2_table = NaN(20,400);  III_16_QLM3_table = NaN(20,400);  

% create auxiliary vectors for results, aux1 for first root, aux2 for
% second root 
aux1_vector = NaN(20,1);   
aux2_vector = NaN(20,1);   

% how close to zero must an eigenvalue be to consider it zero? 
eig_near_zero = 10^-5;

%how close to zero must a derivative be to consider it zero? 
deriv_near_zero = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for r = 1:3         % loop ND (count with r; j is counter for M)
    
    if r == 1         % ND = 1.5
        ND = ND_vector(1);

    else         % N = 8 or 16

        if r == 2         % ND = 8
            ND = ND_vector(2);    

        else         % ND = 16
            ND = ND_vector(3);

        end;         % ND = 8
    end;         % ND = 16
        
    for exp = 2:3        % loop Exp

        % Display Model_StationaryPoint, ND, display experiment number
        disp (['ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m']);

        % Parameter set for currently calculated experiment number
        mu_par = mu_vector(exp-1);
        Phi_P_par = Phi_P_vector(exp-1);
        k_par = k_vector(exp-1);
        g_par = g_vector(exp-1);
        eps_par = eps_vector(exp-1);
        Phi_Z_par = Phi_Z_vector(exp-1);
        beta_par = beta_vector(exp-1);
        m_r_par = m_r_vector(exp-1);
        Omega_par = Omega_vector(exp-1);
        gamma_par = gamma_vector(exp-1);
        Phi_P_opt_par = Phi_P_opt_vector(exp-1);
        Phi_Z_opt_par = Phi_Z_opt_vector(exp-1);

        for j= 1:400         % M: layer depth from 1m to 400 m, calculate terms and variables

            % fill in lines for M, ND, exp
            aux1_vector(posM,1) = j;   aux2_vector(posM,1) = j;   
            aux1_vector(posND,1) = ND;  aux2_vector(posND,1) = ND;
            aux1_vector(posExp,1) = exp; aux2_vector(posExp,1) = exp;
            % calculate auxiliary terms 
                
            % G_1 
            aux1_vector(posG1,1) = (Phi_Z_par*M(j)+m_r_par)/(beta_par*M(j));
            aux2_vector(posG1,1) = (Phi_Z_par*M(j)+m_r_par)/(beta_par*M(j));
            
            %xi_1
            aux1_vector(posxi1,1) = g_par*aux1_vector(posG1,1)/(eps_par*(g_par - eps_par*aux1_vector(posG1,1)));
            aux2_vector(posxi1,1) = g_par*aux2_vector(posG1,1)/(eps_par*(g_par - eps_par*aux2_vector(posG1,1)));
            
            if aux1_vector(posxi1,1) < 0 
                disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  xi1 negative ']);
                countNotifications = countNotifications + 1;
            end;
            if aux2_vector(posxi1,1) < 0 
                disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  xi1 negative ']);
                countNotifications = countNotifications + 1;
            end;
            
            %xi_3
            aux1_vector(posxi3,1) = Phi_P_par/aux1_vector(posG1,1)*sqrt(aux1_vector(posxi1,1))-m_r_par/(M(j)*aux1_vector(posG1,1))*sqrt(aux1_vector(posxi1,1))-Phi_P_opt_par/aux1_vector(posG1,1)*aux1_vector(posxi1,1);
            aux2_vector(posxi3,1) = Phi_P_par/aux2_vector(posG1,1)*sqrt(aux2_vector(posxi1,1))-m_r_par/(M(j)*aux2_vector(posG1,1))*sqrt(aux2_vector(posxi1,1))-Phi_P_opt_par/aux2_vector(posG1,1)*aux2_vector(posxi1,1);

            %xi_4
            aux1_vector(posxi4,1) = -mu_par*sqrt(aux1_vector(posxi1,1))+(1-beta_par)*Omega_par*mu_par*sqrt(aux1_vector(posxi1,1))+Phi_Z_par*Omega_par*mu_par/aux1_vector(posG1,1)*sqrt(aux1_vector(posxi1,1));
            aux2_vector(posxi4,1) = -mu_par*sqrt(aux2_vector(posxi1,1))+(1-beta_par)*Omega_par*mu_par*sqrt(aux2_vector(posxi1,1))+Phi_Z_par*Omega_par*mu_par/aux2_vector(posG1,1)*sqrt(aux2_vector(posxi1,1));

            %xi_5
            aux1_vector(posxi5,1) = gamma_par*Phi_P_par*sqrt(aux1_vector(posxi1,1))-((1-beta_par)*aux1_vector(posG1,1)+Phi_Z_par)*Omega_par*aux1_vector(posxi3,1)+m_r_par/M(j)*aux1_vector(posND,1);
            aux2_vector(posxi5,1) = gamma_par*Phi_P_par*sqrt(aux2_vector(posxi1,1))-((1-beta_par)*aux2_vector(posG1,1)+Phi_Z_par)*Omega_par*aux2_vector(posxi3,1)+m_r_par/M(j)*aux2_vector(posND,1);

            %xi_6
            aux1_vector(posxi6,1) = -aux1_vector(posxi4,1)*M(j)/m_r_par+k_par-aux1_vector(posxi5,1)*M(j)/m_r_par;
            aux2_vector(posxi6,1) = -aux2_vector(posxi4,1)*M(j)/m_r_par+k_par-aux2_vector(posxi5,1)*M(j)/m_r_par;
            
            %xi_7
            aux1_vector(posxi7,1) = -k_par*aux1_vector(posxi5,1)*M(j)/m_r_par;
            aux2_vector(posxi7,1) = -k_par*aux2_vector(posxi5,1)*M(j)/m_r_par;

            %until here aux 1 and aux 2 are the same 
            
            % calculate variables 
            % have MATLAB find roots of N numerically; p = polynom N

            
            p = [1 aux1_vector(posxi6,1) aux1_vector(posxi7,1)];
           
            % MATLAB finds roots via eigen values 
            rootsN_QLM_vector = roots(p);
            
            
           %P
           aux1_vector(posP,1) = sqrt(aux1_vector(posxi1,1));
           aux2_vector(posP,1) = sqrt(aux2_vector(posxi1,1));
           
           if aux1_vector(posxi1,1) < 0 
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  P complex ']);
           end;
           if aux2_vector(posxi1,1) < 0 
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  P complex ']);
           end;
           
           %G
           
           aux1_vector(posG,1) = aux1_vector(posG1,1);
           aux2_vector(posG,1) = aux2_vector(posG1,1);

           if aux1_vector(posG,1) < 0
              disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  G negative ']);
           end;
           if aux1_vector(posG,1) > g_par
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  G > g ']); 
           end;
           
           if aux2_vector(posG,1) < 0
              disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  G negative ']);
           end;
           if aux2_vector(posG,1) > g_par
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  G > g ']); 
           end;
           
           %N_1 for QLM_2
           if abs(imag(rootsN_QLM_vector(1))) > eig_near_zero
               aux1_vector(posN,1) = NaN;
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N: complex root ']); 
           else 
               if real(rootsN_QLM_vector(1)) < 0
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N: negative root ']); 
                   aux1_vector(posN,1) = NaN;
               else
                   aux1_vector(posN,1) = real(rootsN_QLM_vector(1));
                   if aux1_vector(posN,1) > aux1_vector(posND,1) 
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N > ND ']); 
                   end;
               end;
           end;
           
           %N_2 for QLM_3
           if abs(imag(rootsN_QLM_vector(2))) > eig_near_zero
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N: complex root ']); 
               aux2_vector(posN,1) = NaN;
           else
               if real(rootsN_QLM_vector(2)) < 0
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N: negative root ']); 
                   aux2_vector(posN,1) = NaN;
               else
               aux2_vector(posN,1) = real(rootsN_QLM_vector(2));
                   if aux2_vector(posN,1) > aux2_vector(posND,1) 
                   disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N > ND ']); 
                   end;
               end;
           end;
           
           %u
           aux1_vector(posu,1) = aux1_vector(posN,1)/(k_par+aux1_vector(posN,1));
           aux2_vector(posu,1) = aux2_vector(posN,1)/(k_par+aux2_vector(posN,1));
           if aux1_vector(posu,1) < 0 
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u negative ']); 
           end;
           if aux2_vector(posu,1) < 0 
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u negative ']); 
           end;
           if aux1_vector(posu,1) >= 1
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u >= 1 ']); 
           end;
           if aux2_vector(posu,1) >= 1
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u >= 1 ']); 
           end;
           
           %FN
           aux1_vector(posFN,1) = m_r_par/M(j)*(aux1_vector(posND,1)-aux1_vector(posN,1));
           aux2_vector(posFN,1) = m_r_par/M(j)*(aux2_vector(posND,1)-aux2_vector(posN,1));
           if aux1_vector(posFN,1) < 0 
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  ND < 0 ']); 
           end;
           if aux2_vector(posFN,1) < 0 
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  ND < 0 ']); 
           end;
           if aux1_vector > 60
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  ND > 60 ']); 
           end;
           if aux2_vector > 60
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  ND > 60 ']); 
           end;
          
           %Z for QLM_2
           aux1_vector(posZ,1) = mu_par/aux1_vector(posG1,1)*aux1_vector(posu,1)*sqrt(aux1_vector(posxi1,1))-aux1_vector(posxi3,1);
           aux2_vector(posZ,1) = mu_par/aux2_vector(posG1,1)*aux2_vector(posu,1)*sqrt(aux2_vector(posxi1,1))-aux2_vector(posxi3,1);

           if aux1_vector(posZ,1) < 0 
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Z negative ']); 
           end;
           if aux2_vector(posZ,1) < 0 
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Z negative ']); 
           end;
           
           %test derivations
           
           % dotN: aux1_vector(posdotN,1)
           aux1_vector(posdotN,1)= (-mu_par*aux1_vector(posu,1)+gamma_par*Phi_P_par)*aux1_vector(posP,1)+((1-beta_par)*aux1_vector(posG,1)+Phi_Z_par)*Omega_par*aux1_vector(posZ,1)+aux1_vector(posFN,1);
           aux2_vector(posdotN,1)= (-mu_par*aux2_vector(posu,1)+gamma_par*Phi_P_par)*aux2_vector(posP,1)+((1-beta_par)*aux2_vector(posG,1)+Phi_Z_par)*Omega_par*aux2_vector(posZ,1)+aux2_vector(posFN,1);

           % dotN not almost zero
           if abs(aux1_vector(posdotN,1)) > deriv_near_zero
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | N might not be stationary at this depth']);
           end;         % dotN not zero
           if abs(aux2_vector(posdotN,1)) > deriv_near_zero
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | N might not be stationary at this depth']);
           end;         % dotN not zero

           % dotP: aux1_vector(posdotP,1)
           aux1_vector(posdotP,1)= (mu_par*aux1_vector(posu,1)-Phi_P_par-(m_r_par/M(j)))*aux1_vector(posP,1)-aux1_vector(posG,1)*aux1_vector(posZ,1)-Phi_P_opt_par*(aux1_vector(posP,1))^2;
           aux2_vector(posdotP,1)= (mu_par*aux2_vector(posu,1)-Phi_P_par-(m_r_par/M(j)))*aux2_vector(posP,1)-aux2_vector(posG,1)*aux2_vector(posZ,1)-Phi_P_opt_par*(aux2_vector(posP,1))^2;

           % dotP not almost zero
           if abs(aux1_vector(posdotP,1)) > deriv_near_zero
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | P might not be stationary at this depth']);
           end;
           % dotP not almost zero
           if abs(aux2_vector(posdotP,1)) > deriv_near_zero
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | P might not be stationary at this depth']);
           end;

           % dotZ = aux1_vector(posdotZ,1)
           aux1_vector(posdotZ,1)= (beta_par*aux1_vector(posG,1)-Phi_Z_par-(m_r_par/M(j)))*aux1_vector(posZ,1);
           aux2_vector(posdotZ,1)= (beta_par*aux2_vector(posG,1)-Phi_Z_par-(m_r_par/M(j)))*aux2_vector(posZ,1);

           % dotZ not almost zero
           if abs(aux1_vector(posdotZ,1)) > deriv_near_zero
               disp (['QLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | Z might not be stationary at this depth']);
           end;
           if abs(aux2_vector(posdotZ,1)) > deriv_near_zero
               disp (['QLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | Z might not be stationary at this depth']);
           end;
          
           % QLM2 stable? J1 - J9 Jacobi matrix entries 

J1 = -mu_par*(k_par/((k_par+aux1_vector(posN,1))^2))*aux1_vector(posP,1)-m_r_par/M(j);
J2 = (-mu_par*aux1_vector(posu,1)+gamma_par*Phi_P_par)+(1-beta_par)*(2*g_par^2*eps_par*aux1_vector(posP,1)/((g_par+eps_par*aux1_vector(posP,1)^2)^2))*Omega_par*aux1_vector(posZ,1);
J3 = ((1-beta_par)*aux1_vector(posG,1)+Phi_Z_par)*Omega_par; 
J4 = mu_par*(k_par/((k_par+aux1_vector(posN,1))^2))*aux1_vector(posP,1);
J5 = mu_par*aux1_vector(posu,1)-Phi_P_par-m_r_par/M(j)-(2*g_par^2*eps*aux1_vector(posP,1)/((g_par+eps_par*aux1_vector(posP,1)^2)^2))*aux1_vector(posZ,1)-2*Phi_P_opt_par*aux1_vector(posP,1);
J6 = -aux1_vector(posG,1);
J7 = 0;
J8 = beta_par*(2*g_par^2*eps*aux1_vector(posP,1)/((g_par+eps_par*aux1_vector(posP,1)^2)^2))*aux1_vector(posZ,1);
J9 = beta_par*aux1_vector(posG,1)-Phi_Z_par-m_r_par/M(j);

%Jacobi matrix
Jac = [J1 J2 J3; J4 J5 J6; J7 J8 J9];

 % Eigen values of Jacobi matrix: check on NaN entries first
                
                badentry = 0;       % Numbers of NaN entries 
                
                for rowJac = 1:3
                    for columnJac = 1:3
                        entryJac = isnan(Jac(rowJac,columnJac)); 
                        if entryJac == 1 
                            badentry = 1;
                        end; 
                    end;
                end;
                
                if badentry == 1 
                    lambda = [NaN; NaN; NaN];
                else
                    lambda = eig(Jac);
                end;

% test eigenvalues if equilibrium is stable. -1 is default
aux1_vector(posStable,1) = -1;

for lam = 1:3
    if abs(real(lambda(lam))) < eig_near_zero
        aux1_vector(posStable,1) = 0;     % no conclusion possible
    end;
end;

if aux1_vector(posStable,1) == -1
   for lamb = 1:3 
       if real(lambda(lamb)) > 0 
           aux1_vector(posStable,1)= 1;      % unstable 
       end;
   end;
end;

% else: stable

% QLM3 stable? J1 - J9 Jacobi matrix entries 

J1 = -mu_par*(k_par/((k_par+aux2_vector(posN,1))^2))*aux2_vector(posP,1)-m_r_par/M(j);
J2 = (-mu_par*aux2_vector(posu,1)+gamma_par*Phi_P_par)+(1-beta_par)*(2*g_par^2*eps_par*aux2_vector(posP,1)/((g_par+eps_par*aux2_vector(posP,1)^2)^2))*Omega_par*aux2_vector(posZ,1);
J3 = ((1-beta_par)*aux2_vector(posG,1)+Phi_Z_par)*Omega_par; 
J4 = mu_par*(k_par/((k_par+aux2_vector(posN,1))^2))*aux2_vector(posP,1);
J5 = mu_par*aux2_vector(posu,1)-Phi_P_par-m_r_par/M(j)-(2*g_par^2*eps*aux2_vector(posP,1)/((g_par+eps_par*aux2_vector(posP,1)^2)^2))*aux2_vector(posZ,1)-2*Phi_P_opt_par*aux2_vector(posP,1);
J6 = -aux2_vector(posG,1);
J7 = 0;
J8 = beta_par*(2*g_par^2*eps*aux2_vector(posP,1)/((g_par+eps_par*aux2_vector(posP,1)^2)^2))*aux2_vector(posZ,1);
J9 = beta_par*aux2_vector(posG,1)-Phi_Z_par-m_r_par/M(j);

%Jacobi matrix
J = [J1 J2 J3; J4 J5 J6; J7 J8 J9];

 % Eigen values of J: check on NaN entries first
                
                badentry = 0;
              
               
                for rowJac = 1:3
                    for columnJac = 1:3
                        entryJac = isnan(Jac(rowJac,columnJac));
                        if entryJac == 1
                            badentry = 1;
                        end;
                        
                    end;
                end;
                
                if badentry == 1 
                    lambda = [NaN; NaN; NaN];
                else
                    lambda = eig(Jac);
                end;

% test eigenvalues if equilibrium is stable. -1 is default
aux2_vector(posStable,1) = -1;

for lam = 1:3
    if abs(real(lambda(lam))) < eig_near_zero
        aux2_vector(posStable,1) = 0;     % no conclusion possible
    end;
end;

if aux2_vector(posStable,1) == -1
   for lamb = 1:3 
       if real(lambda(lamb)) > 0 
           aux2_vector(posStable,1)= 1;      % unstable 
       end;
   end;
end;

% else stable

           if ND == ND_vector(1)
          
               if exp == 2
                   II_1_QLM2_table(:,j) = aux1_vector;
                   II_1_QLM3_table(:,j) = aux2_vector;
                   
               else         % exp = 3
                   III_1_QLM2_table(:,j) = aux1_vector;
                   III_1_QLM3_table(:,j) = aux2_vector;
               end;         % exp 
           
           else         % ND not 1.5
               
               if ND == ND_vector(2)
                   if exp == 2
                       II_8_QLM2_table(:,j) = aux1_vector;
                       II_8_QLM3_table(:,j) = aux2_vector;
                   else         % exp = 3
                       III_8_QLM2_table(:,j) = aux1_vector;
                       III_8_QLM3_table(:,j) = aux2_vector;

                   end;         % exp 
                   
               else         % ND == ND_vector(3)
                   if exp == 2
                       II_16_QLM2_table(:,j) = aux1_vector;
                       II_16_QLM3_table(:,j) = aux2_vector;

                   else         % exp = 3
                       III_16_QLM2_table(:,j) = aux1_vector;
                       III_16_QLM3_table(:,j) = aux2_vector;

                   end;         % exp 
               end;         % ND not 1.5
           end;         % ND = 1.5
        end;        % j = 1:400
    end;        % loop Exp
end;        % loop ND       
                            

 
%%%%%%%%%%%%%%%%
% QLM plots


figure(1)
subplot(2,2,1) 
hold on;
plot(II_8_QLM2_table(1,:),II_8_QLM2_table(4,:),'b');
plot(II_8_QLM2_table(1,:),II_8_QLM2_table(6,:),'r');
plot(II_8_QLM2_table(1,:),II_8_QLM2_table(8,:),'g');
title('QLM2, Exp.II, N_D=8, N,P,Z');
xlabel('M');
xlim([1,20]);
ylabel('N,P,Z');

subplot(2,2,2)
hold on;
plot(II_16_QLM2_table(1,:),II_16_QLM2_table(4,:),'b');
plot(II_16_QLM2_table(1,:),II_16_QLM2_table(6,:),'r');
plot(II_16_QLM2_table(1,:),II_16_QLM2_table(8,:),'g');
title('QLM2, Exp.II, N_D=16, N,P,Z');
xlabel('M');
xlim([1,400]);
ylabel('N,P,Z');

subplot(2,2,3)
hold on;
plot(II_1_QLM3_table(1,:),II_1_QLM3_table(4,:),'b');
plot(II_1_QLM3_table(1,:),II_1_QLM3_table(6,:),'r');
plot(II_1_QLM3_table(1,:),II_1_QLM3_table(8,:),'g');
title('QLM3, Exp.II, N_D=1, N,P,Z');
xlabel('M');
xlim([1,400]);
ylabel('N,P,Z');
legend('N','P','Z');

subplot(2,2,4)

hold on;
plot(II_8_QLM3_table(1,:),II_8_QLM3_table(4,:),'b');
plot(II_8_QLM3_table(1,:),II_8_QLM3_table(6,:),'r');
plot(II_8_QLM3_table(1,:),II_8_QLM3_table(8,:),'g');
title('QLM3, Exp.II, N_D=8, N,P,Z');
xlabel('M');
xlim([21,400]);
ylabel('N,P,Z');

figure(2)
subplot(2,2,1) 
hold on;
scatter(II_8_QLM2_table(1,:),II_8_QLM2_table(20,:),'r');
title('QLM2, Exp.II, N_D=8, Stability');
xlabel('M');
xlim([1,20]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(2,2,2)
hold on;
scatter(II_16_QLM2_table(1,:),II_16_QLM2_table(20,:),'r');
title('QLM2, Exp.II, N_D=16, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(2,2,3)
hold on;
scatter(II_1_QLM3_table(1,:),II_1_QLM3_table(20,:),'r');
title('QLM3, Exp.II, N_D=1, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(2,2,4)

hold on;
scatter(II_8_QLM3_table(1,:),II_8_QLM3_table(20,:),'r');
title('QLM3, Exp.II, N_D=8, Stability');
xlabel('M');
xlim([21,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

figure(3)
subplot(2,2,1) 
hold on;
plot(III_16_QLM2_table(1,:),III_16_QLM2_table(4,:),'b');
plot(III_16_QLM2_table(1,:),III_16_QLM2_table(6,:),'r');
plot(III_16_QLM2_table(1,:),III_16_QLM2_table(8,:),'g');
title('QLM2, Exp.III, N_D=16, N,P,Z');
xlabel('M');
xlim([65,224]);
ylabel('N,P,Z');

subplot(2,2,2)
hold on;
plot(III_1_QLM3_table(1,:),III_1_QLM3_table(4,:),'b');
plot(III_1_QLM3_table(1,:),III_1_QLM3_table(6,:),'r');
plot(III_1_QLM3_table(1,:),III_1_QLM3_table(8,:),'g');
title('QLM3, Exp.III, N_D=1, N,P,Z');
xlabel('M');
xlim([65,400]);
ylabel('N,P,Z');

subplot(2,2,3)
hold on;
plot(III_8_QLM3_table(1,:),III_8_QLM3_table(4,:),'b');
plot(III_8_QLM3_table(1,:),III_8_QLM3_table(6,:),'r');
plot(III_8_QLM3_table(1,:),III_8_QLM3_table(8,:),'g');
title('QLM3, Exp.III, N_D=8, N,P,Z');
xlabel('M');
xlim([65,400]);
ylabel('N,P,Z');
legend('N','P','Z');

subplot(2,2,4)

hold on;
plot(III_16_QLM3_table(1,:),III_16_QLM3_table(4,:),'b');
plot(III_16_QLM3_table(1,:),III_16_QLM3_table(6,:),'r');
plot(III_16_QLM3_table(1,:),III_16_QLM3_table(8,:),'g');
title('QLM3, Exp.III, N_D=16, N,P,Z');
xlabel('M');
xlim([225,400]);
ylabel('N,P,Z');

figure(4)
subplot(2,2,1) 
scatter(III_16_QLM2_table(1,:),III_16_QLM2_table(20,:),'r');
title('QLM2, Exp.III, N_D=16, Stability');
xlabel('M');
xlim([65,220]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(2,2,2)
hold on;
scatter(III_1_QLM3_table(1,:),III_1_QLM3_table(20,:),'r');
title('QLM2, Exp.II, N_D=16, Stability');
xlabel('M');
xlim([65,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(2,2,3)
hold on;
scatter(III_8_QLM3_table(1,:),III_8_QLM3_table(20,:),'r');
title('QLM3, Exp.II, N_D=1, Stability');
xlabel('M');
xlim([65,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(2,2,4)

hold on;
scatter(III_16_QLM3_table(1,:),II_16_QLM3_table(20,:),'r');
title('QLM3, Exp.III, N_D=16, Stability');
xlabel('M');
xlim([225,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

figure(5)
hold on;
plot(II_8_QLM2_table(1,:),II_8_QLM2_table(4,:),'b');
plot(II_8_QLM2_table(1,:),II_8_QLM2_table(6,:),'r');
plot(II_8_QLM2_table(1,:),II_8_QLM2_table(8,:),'g');
title('QLM2 and 3, Exp.II, N_D=8, N,P,Z');
xlabel('M');
ylabel('N,P,Z');
plot(II_8_QLM3_table(1,:),II_8_QLM3_table(4,:),':b');
plot(II_8_QLM3_table(1,:),II_8_QLM3_table(6,:),'r');
plot(II_8_QLM3_table(1,:),II_8_QLM3_table(8,:),':g');
xlim([1,400]);
legend('N, QLM2','P,QLM2&3','Z,QLM2','N, QLM3','Z,QLM3');



disp 'ready';
disp '********************************';
diary off

               
            
          
   