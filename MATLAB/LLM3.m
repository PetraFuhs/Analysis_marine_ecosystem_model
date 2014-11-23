% LLM3 and LLM4

% Calculating Equilibria NPZ Model
% LLM3 and LLM4 ND = 1.5 or 8 or 16; Parameter set I
% Diploma Thesis P.Fuhs (CAU)
% Matlab Student Version V 7.10.0 (R2010a)


clear all;
clc;
clf;

delete('LLM3_diary'); 
diary on;
diary('LLM3_diary');

disp 'Calculation for LLM3 and LLM4'

% how close to zero must an eigenvalue be to consider it zero? 
eig_near_zero = 10^-6;

%how close to zero must a derivative be to consider it zero? 
deriv_near_zero = 0.0001;

exp = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = NaN(1,400);         % Layer depth from 1 to 400 m, grid = 1m 
for j = 1:400
    M(j)=j;
end; 

% all units as in thesis.
% all domains as in thesis.

% Vector ND values
ND_vector = [1.5, 8, 16];

% Parameter set (experiment I)

mu_par = 2.362;
Phi_P_par = 0.016;
k_par = 0.350;
g_par = 0.743;
eps_par = 3.088;
Phi_Z_par = 0.067;
beta_par = 0.852;
m_r_par = 0.643;
Omega_par = 0.348;
gamma_par = 0.937;

% create helpful result vectors; some of them similar in LLM 3 and 4
% but handled separately for a better overview
% (compare to proof in thesis)

G1_LLM3_vector = NaN(1,400);            G1_LLM4_vector = NaN(1,400);
xi1_LLM3_vector = NaN(1,400);           xi1_LLM4_vector = NaN(1,400);
xi2_LLM3_vector = NaN(1,400);           xi2_LLM4_vector = NaN(1,400);
a1_LLM3_vector = NaN(1,400);            a1_LLM4_vector = NaN(1,400);
b1_LLM3_vector = NaN(1,400);            b1_LLM4_vector = NaN(1,400);
assumption3_LLM3_vector = NaN(1,400);   assumption3_LLM4_vector = NaN(1,400);
N_LLM3_vector = NaN(1,400);             N_LLM4_vector = NaN(1,400);
dotN_LLM3_vector = NaN(1,400);          dotN_LLM4_vector = NaN(1,400);
P_LLM3_vector = NaN(1,400);             P_LLM4_vector = NaN(1,400);
dotP_LLM3_vector = NaN(1,400);          dotP_LLM4_vector = NaN(1,400);
Z_LLM3_vector = NaN(1,400);             Z_LLM4_vector = NaN(1,400);
dotZ_LLM3_vector = NaN(1,400);          dotZ_LLM4_vector = NaN(1,400);
u_LLM3_vector = NaN(1,400);             u_LLM4_vector = NaN(1,400);
G_LLM3_vector = NaN(1,400);             G_LLM4_vector = NaN(1,400);
FN_LLM3_vector = NaN(1,400);            FN_LLM4_vector = NaN(1,400);
stable_LLM3_vector = NaN(1,400);        stable_LLM4_vector = NaN(1,400);

LLM3_countNotifications = 0;
LLM4_countNotifications = 0;

% Create tables with results:
%
% 2 x 3  tables: ND = 1.5,8,16 for LLM3, LLM4
%
% table names: I_NDValue_LLM3_table, e.g.I_8_LLM3_table (same with LLM4)
% 1 in table names refers to ND=1.5
% lines in result matrices:
% line 1: M
% line 2: ND
% line 3: G1
% line 4: xi1
% line 5: xi2
% line 6: a1
% line 7: b1
% line 8: assumption1
% line 9: N
% line 10: dotN
% line 11: P
% line 12: dotP
% line 13: Z
% line 14: dotZ
% line 15: u
% line 16: G
% line 17: FN
% line 18: Stability
%%%%%%%%%%%%%%%%%%%%%%%

I_1_LLM3_table = NaN(18,400);   I_1_LLM4_table = NaN(18,400);
I_8_LLM3_table = NaN(18,400);   I_8_LLM4_table = NaN(18,400);
I_16_LLM3_table = NaN(18,400);  I_16_LLM4_table = NaN(18,400);

% fill in M(j), ND
   
for j = 1:400
    
    I_1_LLM3_table(1,j) = M(j);             I_1_LLM4_table(1,j) = M(j);
    I_1_LLM3_table(2,j) = ND_vector(1);     I_1_LLM4_table(2,j) = ND_vector(1);
    I_8_LLM3_table(1,j) = M(j);             I_8_LLM4_table(1,j) = M(j); 
    I_8_LLM3_table(2,j) = ND_vector(2);     I_8_LLM4_table(2,j) = ND_vector(2);
    I_16_LLM3_table(1,j) = M(j);            I_16_LLM4_table(1,j) = M(j);
    I_16_LLM3_table(2,j) = ND_vector(3);    I_16_LLM4_table(2,j) = ND_vector(3);
    
end; % for j = 1:400 initialize result tables 
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for r = 1:3        %loop ND    
    
    if r == 1 
        ND = ND_vector(1);
         
    else         %  r not 1 

        if r == 2
            ND = ND_vector(2);    

        else         %  r not 2
            ND = ND_vector(3);

        end;         %  r = 2
    end;         %  r = 1
     

    % Display value ND
    disp (['LLM_3, ND: ', num2str(ND)]);
    
    for j = 1:400         %  M: layer depth from 1m to 400 m, calculate terms and variables
        
        
        % calculate auxiliary values and P, G (same for LLM3, LLM4)
        
        
        G1_LLM3_vector(j) = (Phi_Z_par*M(j)+m_r_par)/(beta_par*M(j));
        G1_LLM4_vector(j) = (Phi_Z_par*M(j)+m_r_par)/(beta_par*M(j));

        
        assumption3_LLM3_vector(j) = g_par*beta_par*M(j)-Phi_Z_par*M(j)-m_r_par;
        assumption3_LLM4_vector(j) = g_par*beta_par*M(j)-Phi_Z_par*M(j)-m_r_par;

        if assumption3_LLM3_vector(j) <= 0
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 3 violated ']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;
        if assumption3_LLM4_vector(j) <= 0
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 3 violated ']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;
        
        xi1_LLM3_vector(j) = g_par*G1_LLM3_vector(j)/(eps_par*(g_par-G1_LLM3_vector(j)));
        xi1_LLM4_vector(j) = g_par*G1_LLM4_vector(j)/(eps_par*(g_par-G1_LLM4_vector(j)));
        
        if xi1_LLM3_vector(j) < 0 % comes with violation of assumption 3
            xi1_LLM3_vector(j) = NaN;
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  xi1 negative, set to NaN ']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;   
        else 
            P_LLM3_vector(j) = sqrt(xi1_LLM3_vector(j)); 
        end;
        if xi1_LLM4_vector(j) < 0 % comes with violation of assumption 3
            xi1_LLM4_vector(j) = NaN;
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  xi1 negative, set to NaN ']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;   
        else 
            P_LLM4_vector(j) = sqrt(xi1_LLM4_vector(j)); 
        end;
       

        G_LLM3_vector(j) = G1_LLM3_vector(j);
        G_LLM4_vector(j) = G1_LLM4_vector(j);

               
        xi2_LLM3_vector(j) = (((1-beta_par)*G1_LLM3_vector(j)+Phi_Z_par)*Omega_par*sqrt(xi1_LLM3_vector(j))/G1_LLM3_vector(j));     %xi1 is greater zero due to ass.1
        xi2_LLM4_vector(j) = (((1-beta_par)*G1_LLM4_vector(j)+Phi_Z_par)*Omega_par*sqrt(xi1_LLM4_vector(j))/G1_LLM4_vector(j));     %xi1 is greater zero due to ass.1

        a1_LLM3_vector(j)= ((M(j)/m_r_par)*(sqrt(xi1_LLM3_vector(j))*(mu_par-gamma_par*Phi_P_par)-xi2_LLM3_vector(j)*(mu_par-Phi_P_par-m_r_par/M(j))+m_r_par/M(j)*(k_par-ND)));
        a1_LLM4_vector(j)= ((M(j)/m_r_par)*(sqrt(xi1_LLM4_vector(j))*(mu_par-gamma_par*Phi_P_par)-xi2_LLM4_vector(j)*(mu_par-Phi_P_par-m_r_par/M(j))+m_r_par/M(j)*(k_par-ND)));

        b1_LLM3_vector(j)= (k_par*M(j)/m_r_par*(-gamma_par*Phi_P_par*sqrt(xi1_LLM3_vector(j))+xi2_LLM3_vector(j)*(Phi_P_par+m_r_par/M(j)))-k_par*ND);
        b1_LLM4_vector(j)= (k_par*M(j)/m_r_par*(-gamma_par*Phi_P_par*sqrt(xi1_LLM4_vector(j))+xi2_LLM4_vector(j)*(Phi_P_par+m_r_par/M(j)))-k_par*ND);

        if b1_LLM3_vector(j) >= 1/4*a1_LLM3_vector(j) 
            LLM3_countNotifications = LLM3_countNotifications +1;
            N_LLM3_vector(j) = NaN;
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 4 violated ']);            
        end;
        if b1_LLM4_vector(j) >= 1/4*a1_LLM4_vector(j) 
            LLM4_countNotifications = LLM4_countNotifications +1;
            N_LLM4_vector(j) = NaN;
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 4 violated ']);            
        end;
        
        if a1_LLM3_vector(j) >= 0 
            if b1_LLM3_vector(j) >= 0
                LLM3_countNotifications = LLM3_countNotifications + 1;
                N_LLM3_vector(j) = NaN;
                disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 5 violated, ']);            
            end;
        end;
        if a1_LLM4_vector(j) >= 0 
            if b1_LLM4_vector(j) >= 0
                LLM4_countNotifications = LLM4_countNotifications + 1;
                N_LLM4_vector(j) = NaN;
                disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 5 violated, ']);            
            end;
        end;
        
        
        % N, u, FN, Z, dotN, dotP, dotZ for LLM3
        N_LLM3_vector(j)= -0.5*a1_LLM3_vector(j)+0.5*sqrt((a1_LLM3_vector(j))^2-4*b1_LLM3_vector(j));
        
        if N_LLM3_vector(j) < 0 
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N is negative ']);           
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;
        
        if N_LLM3_vector(j) > ND
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |   N exceeds ND ']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;
       
        u_LLM3_vector(j) = N_LLM3_vector(j)/(N_LLM3_vector(j)+k_par);
        
        if u_LLM3_vector(j) < 0 
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u negative ']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;
        
        if u_LLM3_vector(j) >= 1
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u >= 1 ']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;
        
        FN_LLM3_vector(j) = m_r_par/M(j)*(ND-N_LLM3_vector(j));
        
        if FN_LLM3_vector(j) < 0 
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  FN negative ']);            
        end;
        
        if FN_LLM3_vector(j) > 60
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  FN too big ']);            
        end;
        
        Z_LLM3_vector(j) = sqrt(xi1_LLM3_vector(j))/G1_LLM3_vector(j)*(mu_par*u_LLM3_vector(j)-Phi_P_par-m_r_par/M(j));
        
        if Z_LLM3_vector(j) < 0 
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Z is negative']);            
           LLM3_countNotifications = LLM3_countNotifications + 1;
        end;
        
        dotN_LLM3_vector(j)= (-mu_par*u_LLM3_vector(j)+gamma_par*Phi_P_par)*P_LLM3_vector(j)+((1-beta_par)*G_LLM3_vector(j)+Phi_Z_par)*Omega_par*Z_LLM3_vector(j)+FN_LLM3_vector(j);

        if  abs(dotN_LLM3_vector(j)) > deriv_near_zero
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N might not be stationary at this depth']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;        % N not stationary

        dotP_LLM3_vector(j)= (mu_par*u_LLM3_vector(j)-Phi_P_par-m_r_par/M(j))*P_LLM3_vector(j)-G_LLM3_vector(j)*Z_LLM3_vector(j);

        if  abs(dotP_LLM3_vector(j)) > deriv_near_zero
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | Notification:P might not be stationary at this depth']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;    % P not stationary

        dotZ_LLM3_vector(j) = (beta_par*G_LLM3_vector(j)-Phi_Z_par-m_r_par/M(j))*Z_LLM3_vector(j);

        if  abs(dotZ_LLM3_vector(j)) > deriv_near_zero
            disp (['LLM3, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Z might not be stationary at this depth']);            
            LLM3_countNotifications = LLM3_countNotifications + 1;
        end;     % Z not stationary
        
        % stability 
                % stable? J1 - J9 Jacobian matrix entries 

                J1 = -mu_par*(k_par/((k_par+N_LLM3_vector(j))^2))*P_LLM3_vector(j)-m_r_par/M(j);
                J2 = (-mu_par*u_LLM3_vector(j)+gamma_par*Phi_P_par)+(1-beta_par)*((2*g_par^2*eps_par*P_LLM3_vector(j))/((g_par+eps_par*(P_LLM3_vector(j))^2)^2))*Omega_par*Z_LLM3_vector(j);
                J3 = ((1-beta_par)*G_LLM3_vector(j)+Phi_Z_par)*Omega_par; 
                J4 = mu_par*(k_par/((k_par+N_LLM3_vector(j))^2))*P_LLM3_vector(j);
                J5 = mu_par*u_LLM3_vector(j)-Phi_P_par-m_r_par/M(j)-((2*g_par^2*eps*P_LLM3_vector(j))/((g_par+eps_par*(P_LLM3_vector(j))^2)^2))*Z_LLM3_vector(j);
                J6 = -G_LLM3_vector(j);
                J7 = 0;
                J8 = beta_par*((2*g_par^2*eps*P_LLM3_vector(j))/(((g_par+eps_par*(P_LLM3_vector(j))^2)^2)))*Z_LLM3_vector(j);
                J9 = beta_par*G_LLM3_vector(j)-Phi_Z_par-m_r_par/M(j);

                %Jacobian matrix
                Jac = [J1 J2 J3; J4 J5 J6; J7 J8 J9];

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

                % test eigenvalues if equilibrium is stable. 2 is default 
                % value for "not yet tested" 
                
                if badentry == 0
                stable_LLM3_vector(j) = 2;

                    for lam = 1:3
                        if abs(real(lambda(lam))) < eig_near_zero
                            stable_LLM3_vector(j) = 0;     % no conclusion possible
                        end;
                    end;
                
                    if stable_LLM3_vector(j) == 2
                       for lamb = 1:3 
                           if real(lambda(lamb)) > 0 
                               stable_LLM3_vector(j)= 1;      % unstable 
                           end;
                       end;
                    end;

                    if stable_LLM3_vector(j) == 2
                       stable_LLM3_vector(j) = -1;             % stable
                    end;
                end;

        % calculate N, u, FN, Z, dotN, dotP, dotZ for LLM4 
        
       
        N_LLM4_vector(j)= -0.5*a1_LLM3_vector(j)-0.5*sqrt((a1_LLM3_vector(j))^2-4*b1_LLM3_vector(j));
        
        if N_LLM4_vector(j) < 0 
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N is negative']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;
        
        u_LLM4_vector(j) = N_LLM4_vector(j)/(N_LLM4_vector(j)+k_par);
        
        if u_LLM4_vector(j) < 0 
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u negative']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;
        
        if u_LLM4_vector(j) >= 1
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u >= 1']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;
        
        FN_LLM4_vector(j) = m_r_par/M(j)*(ND-N_LLM4_vector(j));
        
        if FN_LLM4_vector(j) < 0
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | FN negative']);            
        end;
        
        if FN_LLM4_vector(j) > 60
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | FN too big']);            
        end;
        
         Z_LLM4_vector(j) = beta_par*M(j)/(Phi_Z_par*M(j)+m_r_par)*(mu_par*u_LLM4_vector(j)-Phi_P_par-m_r_par/M(j))*sqrt(g_par*(Phi_Z_par*M(j)+m_r_par)/(eps_par*(g_par*beta_par*M(j)-Phi_Z_par*M(j)-m_r_par)));
        
        if Z_LLM4_vector(j) < 0 
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | Z is negative']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;
        
        dotN_LLM4_vector(j)= (-mu_par*u_LLM4_vector(j)+gamma_par*Phi_P_par)*P_LLM4_vector(j)+((1-beta_par)*G_LLM4_vector(j)+Phi_Z_par)*Omega_par*Z_LLM4_vector(j)+FN_LLM4_vector(j);

        if  abs(dotN_LLM4_vector(j)) > deriv_near_zero
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | N might be not stationary at this depth']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;        % N not stationary

        dotP_LLM4_vector(j)= (mu_par*u_LLM4_vector(j)-Phi_P_par-m_r_par/M(j))*P_LLM4_vector(j)-G_LLM4_vector(j)*Z_LLM4_vector(j);

        if  abs(dotP_LLM4_vector(j)) > deriv_near_zero
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | P might be not stationary at this depth']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;    % P not stationary

        dotZ_LLM4_vector(j) = (beta_par*G_LLM4_vector(j)-Phi_Z_par-m_r_par/M(j))*Z_LLM4_vector(j);

        if  abs(dotZ_LLM4_vector(j)) > deriv_near_zero
            disp (['LLM4, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Z might not be stationary at this depth']);            
            LLM4_countNotifications = LLM4_countNotifications + 1;
        end;     % Z not stationary

        % stability 
                % stable? J1 - J9 Jacobian matrix entries 

                J1 = -mu_par*(k_par/((k_par+N_LLM3_vector(j))^2))*P_LLM3_vector(j)-m_r_par/M(j);
                J2 = (-mu_par*u_LLM3_vector(j)+gamma_par*Phi_P_par)+(1-beta_par)*((2*g_par^2*eps_par*P_LLM3_vector(j))/((g_par+eps_par*(P_LLM3_vector(j))^2)^2))*Omega_par*Z_LLM3_vector(j);
                J3 = ((1-beta_par)*G_LLM3_vector(j)+Phi_Z_par)*Omega_par; 
                J4 = mu_par*(k_par/((k_par+N_LLM3_vector(j))^2))*P_LLM3_vector(j);
                J5 = mu_par*u_LLM3_vector(j)-Phi_P_par-m_r_par/M(j)-((2*g_par^2*eps*P_LLM3_vector(j))/((g_par+eps_par*(P_LLM3_vector(j))^2)^2))*Z_LLM3_vector(j);
                J6 = -G_LLM3_vector(j);
                J7 = 0;
                J8 = beta_par*((2*g_par^2*eps*P_LLM3_vector(j))/(((g_par+eps_par*(P_LLM3_vector(j))^2)^2)))*Z_LLM3_vector(j);
                J9 = beta_par*G_LLM3_vector(j)-Phi_Z_par-m_r_par/M(j);


                %Jacobian matrix
                Jac = [J1 J2 J3; J4 J5 J6; J7 J8 J9];

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

                % test eigenvalues if equilibrium is stable. 2 is default 
                % value for "not yet tested" 
                
                
                if badentry == 0
                    stable_LLM4_vector(j) = 2;
                    for lam = 1:3
                        if abs(real(lambda(lam))) < eig_near_zero
                            stable_LLM4_vector(j) = 0;     % no conclusion possible
                        end;
                    end;
                
                    if stable_LLM4_vector(j) == 2
                       for lamb = 1:3 
                           if real(lambda(lamb)) > 0 
                               stable_LLM4_vector(j)= 1;      % unstable 
                           end;
                       end;
                    end;

                    if stable_LLM4_vector(j) == 2
                       stable_LLM4_vector(j) = -1;             % stable
                    end;
                end;
     if ND == ND_vector(1)

                    I_1_LLM3_table(3,j) = G1_LLM3_vector(j);
                    I_1_LLM3_table(4,j) = xi1_LLM3_vector(j);
                    I_1_LLM3_table(5,j) = xi2_LLM3_vector(j);
                    I_1_LLM3_table(6,j) = a1_LLM3_vector(j);
                    I_1_LLM3_table(7,j) = b1_LLM3_vector(j);
                    I_1_LLM3_table(8,j) = assumption3_LLM3_vector(j);                        
                    I_1_LLM3_table(9,j) = N_LLM3_vector(j);
                    I_1_LLM3_table(10,j) = dotN_LLM3_vector(j);
                    I_1_LLM3_table(11,j) = P_LLM3_vector(j);
                    I_1_LLM3_table(12,j) = dotP_LLM3_vector(j);
                    I_1_LLM3_table(13,j) = Z_LLM3_vector(j);
                    I_1_LLM3_table(14,j) = dotZ_LLM3_vector(j);
                    I_1_LLM3_table(15,j) = u_LLM3_vector(j);
                    I_1_LLM3_table(16,j) = G_LLM3_vector(j);
                    I_1_LLM3_table(17,j) = FN_LLM3_vector(j);
                    I_1_LLM3_table(18,j) = stable_LLM3_vector(j);
                    
                    I_1_LLM4_table(3,j) = G1_LLM4_vector(j);
                    I_1_LLM4_table(4,j) = xi1_LLM4_vector(j);
                    I_1_LLM4_table(5,j) = xi2_LLM4_vector(j);
                    I_1_LLM4_table(6,j) = a1_LLM4_vector(j);
                    I_1_LLM4_table(7,j) = b1_LLM4_vector(j);
                    I_1_LLM4_table(8,j) = assumption3_LLM4_vector(j);                        
                    I_1_LLM4_table(9,j) = N_LLM4_vector(j);
                    I_1_LLM4_table(10,j) = dotN_LLM4_vector(j);
                    I_1_LLM4_table(11,j) = P_LLM4_vector(j);
                    I_1_LLM4_table(12,j) = dotP_LLM4_vector(j);
                    I_1_LLM4_table(13,j) = Z_LLM4_vector(j);
                    I_1_LLM4_table(14,j) = dotZ_LLM4_vector(j);
                    I_1_LLM4_table(15,j) = u_LLM4_vector(j);
                    I_1_LLM4_table(16,j) = G_LLM4_vector(j);
                    I_1_LLM4_table(17,j) = FN_LLM4_vector(j);
                    I_1_LLM3_table(18,j) = stable_LLM4_vector(j);


       else
           if ND == ND_vector(2)

                    I_8_LLM3_table(3,j) = G1_LLM3_vector(j);
                    I_8_LLM3_table(4,j) = xi1_LLM3_vector(j);
                    I_8_LLM3_table(5,j) = xi2_LLM3_vector(j);
                    I_8_LLM3_table(6,j) = a1_LLM3_vector(j);
                    I_8_LLM3_table(7,j) = b1_LLM3_vector(j);
                    I_8_LLM3_table(8,j) = assumption3_LLM3_vector(j);                        
                    I_8_LLM3_table(9,j) = N_LLM3_vector(j);
                    I_8_LLM3_table(10,j) = dotN_LLM3_vector(j);
                    I_8_LLM3_table(11,j) = P_LLM3_vector(j);
                    I_8_LLM3_table(12,j) = dotP_LLM3_vector(j);
                    I_8_LLM3_table(13,j) = Z_LLM3_vector(j);
                    I_8_LLM3_table(14,j) = dotZ_LLM3_vector(j);
                    I_8_LLM3_table(15,j) = u_LLM3_vector(j);
                    I_8_LLM3_table(16,j) = G_LLM3_vector(j);
                    I_8_LLM3_table(17,j) = FN_LLM3_vector(j);
                    I_8_LLM3_table(18,j) = stable_LLM3_vector(j);

                    I_8_LLM4_table(3,j) = G1_LLM4_vector(j);
                    I_8_LLM4_table(4,j) = xi1_LLM4_vector(j);
                    I_8_LLM4_table(5,j) = xi2_LLM4_vector(j);
                    I_8_LLM4_table(6,j) = a1_LLM4_vector(j);
                    I_8_LLM4_table(7,j) = b1_LLM4_vector(j);
                    I_8_LLM4_table(8,j) = assumption3_LLM4_vector(j);                        
                    I_8_LLM4_table(9,j) = N_LLM4_vector(j);
                    I_8_LLM4_table(10,j) = dotN_LLM4_vector(j);
                    I_8_LLM4_table(11,j) = P_LLM4_vector(j);
                    I_8_LLM4_table(12,j) = dotP_LLM4_vector(j);
                    I_8_LLM4_table(13,j) = Z_LLM4_vector(j);
                    I_8_LLM4_table(14,j) = dotZ_LLM4_vector(j);
                    I_8_LLM4_table(15,j) = u_LLM4_vector(j);
                    I_8_LLM4_table(16,j) = G_LLM4_vector(j);
                    I_8_LLM4_table(17,j) = FN_LLM4_vector(j);
                    I_8_LLM4_table(18,j) = stable_LLM4_vector(j);

           else         %  ND==3                       

                    I_1_LLM3_table(3,j) = G1_LLM3_vector(j);
                    I_16_LLM3_table(4,j) = xi1_LLM3_vector(j);
                    I_16_LLM3_table(5,j) = xi2_LLM3_vector(j);
                    I_16_LLM3_table(6,j) = a1_LLM3_vector(j);
                    I_16_LLM3_table(7,j) = b1_LLM3_vector(j);
                    I_16_LLM3_table(8,j) = assumption3_LLM3_vector(j);                        
                    I_16_LLM3_table(9,j) = N_LLM3_vector(j);
                    I_16_LLM3_table(10,j) = dotN_LLM3_vector(j);
                    I_16_LLM3_table(11,j) = P_LLM3_vector(j);
                    I_16_LLM3_table(12,j) = dotP_LLM3_vector(j);
                    I_16_LLM3_table(13,j) = Z_LLM3_vector(j);
                    I_16_LLM3_table(14,j) = dotZ_LLM3_vector(j);
                    I_16_LLM3_table(15,j) = u_LLM3_vector(j);
                    I_16_LLM3_table(16,j) = G_LLM3_vector(j);
                    I_16_LLM3_table(17,j) = FN_LLM3_vector(j);
                    I_16_LLM3_table(18,j) = stable_LLM3_vector(j);

                    
                    I_1_LLM4_table(3,j) = G1_LLM4_vector(j);
                    I_16_LLM4_table(4,j) = xi1_LLM4_vector(j);
                    I_16_LLM4_table(5,j) = xi2_LLM4_vector(j);
                    I_16_LLM4_table(6,j) = a1_LLM4_vector(j);
                    I_16_LLM4_table(7,j) = b1_LLM4_vector(j);
                    I_16_LLM4_table(8,j) = assumption3_LLM4_vector(j);                        
                    I_16_LLM4_table(9,j) = N_LLM4_vector(j);
                    I_16_LLM4_table(10,j) = dotN_LLM4_vector(j);
                    I_16_LLM4_table(11,j) = P_LLM4_vector(j);
                    I_16_LLM4_table(12,j) = dotP_LLM4_vector(j);
                    I_16_LLM4_table(13,j) = Z_LLM4_vector(j);
                    I_16_LLM4_table(14,j) = dotZ_LLM4_vector(j);
                    I_16_LLM4_table(15,j) = u_LLM4_vector(j);
                    I_16_LLM4_table(16,j) = G_LLM4_vector(j);
                    I_16_LLM4_table(17,j) = FN_LLM4_vector(j);
                    I_16_LLM4_table(18,j) = stable_LLM4_vector(j);

           end;         %  if ND=2
      end;         %  if ND=1
      
    end;     % j = 1:400      % loop ND
end;

disp ('Notifications regarding LLM3: ')
disp (LLM3_countNotifications);
disp ('Notifications regarding LLM4: ')
disp (LLM4_countNotifications);
disp 'ready';
diary off;

figure(1)

subplot(1,2,1)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(8,:)); 
title('LLM3, Assumption 3');
refline(0,0);
xlabel('M');
xlim([1,400]);
ylabel('Assumption 3');

subplot(1,2,2)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(8,:)); 
scatter(I_1_LLM3_table(1,:),I_1_LLM3_table(8,:),'r'); 
title('LLM3, Assumption 3, Detail');
xlim([1,3]);
refline(0,0);
ylim([-0.2,1.1]);
xlabel('M');
ylabel('Assumption 3');

figure(2)

subplot(3,2,1)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(9,:),'b'); 
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(9,:),'r'); 
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(9,:),'g'); 
title('LLM3, N');
xlabel('M');
xlim([1,400]);
ylabel('N');
legend('N_D = 1.5','N_D = 8','N_D = 16');

subplot(3,2,2)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(9,:),'b'); 
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(9,:),'r'); 
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(9,:),'g'); 
title('LLM3, N, Detail');
xlim([0,10]);
ylim([0,12.1]);
xlabel('M');
ylabel('N');

subplot(3,2,3)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(11,:),'b'); 
title('LLM3, P');
xlabel('M');
xlim([1,400]);
ylabel('P');

subplot(3,2,4)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(11,:),'b'); 
title('LLM3, P, Detail');
xlim([0,50]);
xlabel('M');
ylabel('P');

subplot(3,2,5)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(13,:),'b'); 
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(13,:),'r'); 
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(13,:),'g'); 
title('LLM3, Z');
xlabel('M');
xlim([1,400]);
ylabel('Z');

subplot(3,2,6)
hold on;
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(13,:),'b'); 
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(13,:),'r'); 
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(13,:),'g'); 
title('LLM3, Z, Detail');
xlim([1,100]);
xlabel('M');
ylabel('Z');

figure(3)
subplot(3,2,1)
hold on; 
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(9,:),'b');
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(11,:),'r');
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(13,:),'g');
title('LLM3, ND=1.5, N,P,Z on M');
xlabel('M');
xlim([1,400]);
ylabel('N,P,Z');

subplot(3,2,2)
hold on; 
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(9,:),'b');
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(11,:),'r');
plot(I_1_LLM3_table(1,:),I_1_LLM3_table(13,:),'g');
title('LLM3, ND=1.5, N,P,Z, Detail');
xlabel('M');
xlim([1,30]);
ylabel('N,P,Z');
ylim([0,0.75]);

subplot(3,2,3)
hold on; 
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(9,:),'b');
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(11,:),'r');
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(13,:),'g');
title('LLM3, ND=8, N,P,Z');
xlabel('M');
xlim([1,400]);
ylabel('N,P,Z');

subplot(3,2,4)
hold on; 
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(9,:),'b');
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(11,:),'r');
plot(I_8_LLM3_table(1,:),I_8_LLM3_table(13,:),'g');
title('LLM3, ND=8, N,P,Z on M, Detail');
xlabel('M');
xlim([1,30]);
ylabel('N,P,Z');
ylim([0,5]);

subplot(3,2,5)
hold on; 
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(9,:),'b');
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(11,:),'r');
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(13,:),'g');
title('LLM3, ND=16, N,P,Z on M');
xlabel('M');
xlim([1,400]);
ylabel('N,P,Z');
legend('N','P', 'Z');

subplot(3,2,6)
hold on; 
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(9,:),'b');
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(11,:),'r');
plot(I_16_LLM3_table(1,:),I_16_LLM3_table(13,:),'g');
title('LLM3, ND=16, N,P,Z on M, Detail');
xlabel('M');
xlim([1,30]);
ylabel('N,P,Z')
ylim([0,13]);

figure(4)

subplot(1,3,1)
scatter(I_1_LLM3_table(1,:),I_1_LLM3_table(18,:),'r');
title('LLM3, Exp.1, ND=1.5, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(1,3,2)
scatter(I_8_LLM3_table(1,:),I_8_LLM3_table(18,:),'r');
title('LLM3, Exp.1, ND=8, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(1,3,3)
scatter(I_16_LLM3_table(1,:),I_16_LLM3_table(18,:),'r');
title('LLM3, Exp.1, ND=16, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

% subplot(2,3,4)
% scatter(I_1_LLM4_table(1,:),I_1_LLM4_table(18,:));
% title('LLM4, Exp.1, ND=1.5, Stability');
% xlabel('M');
% ylim([-2,2]);
% ylabel('Stability');

% subplot(2,3,5)
% scatter(I_8_LLM4_table(1,:),I_8_LLM4_table(18,:),'r');
% title('LLM4, Exp.1, ND=8, Stability');
% xlabel('M');
% ylim([-2,2]);
% ylabel('Stability');
% 
% subplot(2,3,6)
% scatter(I_16_LLM4_table(1,:),I_16_LLM4_table(18,:),'r');
% title('LLM4, Exp.1, ND=16, Stability');
% xlabel('M');
% ylim([-2,2]);
% ylabel('Stability');
