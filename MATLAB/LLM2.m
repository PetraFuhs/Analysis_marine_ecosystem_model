% LLM2

% Calculating Equilibrium LLM 2 NPZ Model 
% LLM2 ND = 1.5 or 8 or 16; Parameter set I
% Diploma Thesis P.Fuhs (CAU)
% Matlab Student Version V 7.10.0 (R2010a)

clear all;
clc;
clf;

delete('LLM2_diary');
diary on;
diary('LLM2_diary');
countNotifications = 0;

disp 'Calculation for LLM2'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = NaN(1,400);         % M: Layer depth from 1 to 400, grid = 1 
for j = 1:400
    M(j)=j;
end; 

% All units as in thesis.
% All domains as in thesis.
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

% create helpful result vectors

N_LLM2_vector = NaN(1,400);
dotN_LLM2_vector = NaN(1,400);
P_LLM2_vector = NaN(1,400);
dotP_LLM2_vector = NaN(1,400);
Z_LLM2_vector = NaN(1,400);
dotZ_LLM2_vector = NaN(1,400);
assumption1_LLM2_vector = NaN(1,400);
assumption2_LLM2_vector = NaN(1,400);
u_LLM2_vector = NaN(1,400);
G_LLM2_vector = NaN(1,400);
FN_LLM2_vector = NaN(1,400);
stable_LLM2_vector = NaN(1,400);


% Create tables with results:
%
% 3  tables: ND = 1.5,8,16
%
% table names: ExperimentNumber_NDValue_LLM2_table, e.g.I_8_LLM2_table
% ND=1.5 refers to 1 in table names
% lines in result matrices:
% line 1: M
% line 2: ND
% line 3: N
% line 4: dotN
% line 5: P
% line 6: dotP
% line 7: Z
% line 8: dotZ
% line 9: Assumption 1
% line 10: Assumption 2
% line 11: u
% line 12: G
% line 13: FN
% line 14: stability
%

I_1_LLM2_table = NaN(14,400); 
I_8_LLM2_table = NaN(14,400);
I_16_LLM2_table = NaN(14,400);

% how close to zero must an eigenvalue be to consider it zero? 
eig_near_zero = 0.001;

%how close to zero must a derivative be to consider it zero? 
deriv_near_zero = 0.0001;

exp = 1;
% fill in M(j), ND
   
for j = 1:400
    
    I_1_LLM2_table(1,j) = M(j);
    I_1_LLM2_table(2,j) = ND_vector(1);
    I_8_LLM2_table(1,j) = M(j); 
    I_8_LLM2_table(2,j) = ND_vector(2);
    I_16_LLM2_table(1,j) = M(j);
    I_16_LLM2_table(2,j) = ND_vector(3);
    
end; % for j = 1:400 initialize result tables 
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for r = 1:3         % loop ND 
    
    if r == 1 
        ND = ND_vector(1);
         
    else         %  r not 1 

        if r == 2
            ND = ND_vector(2);    

        else         %  r not 2
            ND = ND_vector(3);

        end;         %  r = 2
    end;         %  r = 1
    

    %       Display value ND
    disp (['LLM_2, ND: ', num2str(ND)]);

    for j = 1:400         %  M: layer depth from 1m to 400, calculate terms and variables

    assumption1_LLM2_vector(j) =  mu_par * M(j) - Phi_P_par * M(j) - m_r_par;

    if j == 1         %  avoid zero crossing assumption 1 

        if assumption1_LLM2_vector(j) == 0         %  avoid zero in M=1
        disp (['ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 1 is zero in digit 1']);
        countNotifications = countNotifications + 1;
        keyboard;
        end;        % avoid zero assumption 1 in M = 1

    else         %  avoid zero crossing assumption 1 in M > 1
        if assumption1_LLM2_vector(j)*assumption1_LLM2_vector(j-1) <= 0
        disp (['ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 1 is zero between ', num2str(j-1), ' and ', num2str(j)']);
        countNotifications = countNotifications + 1;
        keyboard;
        end,

    end;         %  avoid zero-crossing assumption 1 

    assumption2_LLM2_vector(j) = gamma_par * Phi_P_par * M(j) - Phi_P_par * M(j) - m_r_par;

        if j == 1         %  avoid zero crossing assumption 2

            if assumption2_LLM2_vector(j) == 0         %  avoid zero in assumption 2 in M = 1
            disp(['Experiment ', num2str(exp), ' Assumption 2 is zero in digit 1']);
            disp (['ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Assumption 2 is zero in digit 1']);
            countNotifications = countNotifications + 1;
            keyboard;
            end;

        else         %  avoid zero crossing assumption 2 in M > 1

            if assumption2_LLM2_vector(j)*assumption2_LLM2_vector(j-1)<= 0
            disp(['LLM2, Experiment ', num2str(exp); 'ND: ', num2str(exp); '  Assumption 2 is zero beween digits ', num2str(j-1), ' and ', num2str(j)]); 
            countNotifications = countNotifications + 1;
            keyboard;
            end,

        end;         %  avoid zero crossing assumption 2

        N_LLM2_vector(j) = k_par * ((Phi_P_par*(M(j))+m_r_par)/(mu_par*M(j)-Phi_P_par*M(j)-m_r_par));
        
        % N > ND ? 
        if N_LLM2_vector(j) > ND
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N exceeds ND ']);
            countNotifications = countNotifications + 1;
        end;
        
        P_LLM2_vector(j) = (m_r_par/(gamma_par*Phi_P_par*M(j)-Phi_P_par*M(j)-m_r_par))*(-ND + k_par*(Phi_P_par*M(j)+m_r_par)/(mu_par*M(j)-Phi_P_par*M(j)-m_r_par));

        Z_LLM2_vector(j) = 0;

        % calculate terms u, FN, G
        % Check: Are derivatives almost zero? 

        u_LLM2_vector(j) = N_LLM2_vector(j)/(N_LLM2_vector(j)+k_par);
 
        % u < 0 ? 
        if u_LLM2_vector(j) < 0 
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u negative ']);
            countNotifications = countNotifications + 1;
        end;

        % u >= 1 ? 
        if u_LLM2_vector(j) >= 1
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u >= 1']);
            countNotifications = countNotifications + 1;
        end;

        G_LLM2_vector(j) = (g_par*eps_par*((P_LLM2_vector(j))^2))/(g_par+eps_par*((P_LLM2_vector(j))^2));

        % G < 0 ? 
        if G_LLM2_vector(j) < 0
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  G negative']);
            countNotifications = countNotifications + 1;
        end;

        % G >= g? 
        if G_LLM2_vector(j) >= g_par
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  G >= g ']);
            countNotifications = countNotifications + 1;
        end;

        
        FN_LLM2_vector(j) = m_r_par/M(j)*(ND-N_LLM2_vector(j));

        % FN < 0 ?
        if FN_LLM2_vector(j) < 0
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  F_N negative']);
            countNotifications = countNotifications + 1;
        end;
        
        % FN > 60 ? 
            if FN_LLM2_vector(j)> 60
                disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  F_N > 60']);
                countNotifications = countNotifications + 1;
            end;   

        dotN_LLM2_vector(j) = (-mu_par*u_LLM2_vector(j)+gamma_par*Phi_P_par)*P_LLM2_vector(j)+((1-beta_par)*G_LLM2_vector(j)+Phi_Z_par)*Omega_par*Z_LLM2_vector(j)+FN_LLM2_vector(j);

        if abs(dotN_LLM2_vector(j)) > deriv_near_zero
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  N might be not stationary at this depth']);
            countNotifications = countNotifications + 1;
        end;

        dotP_LLM2_vector(j) = (mu_par*u_LLM2_vector(j)-Phi_P_par-m_r_par/M(j))*P_LLM2_vector(j)-G_LLM2_vector(j)*Z_LLM2_vector(j);

        if abs(dotP_LLM2_vector(j)) > deriv_near_zero
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  P might not be stationary at this depth']);
            countNotifications = countNotifications + 1;
        end;

        dotZ_LLM2_vector(j) = (beta_par*G_LLM2_vector(j)-Phi_Z_par-m_r_par/M(j))*Z_LLM2_vector(j);

        if abs(dotZ_LLM2_vector(j)) > deriv_near_zero
            disp (['LLM2, ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Z might not be stationary at this depth ']);
            countNotifications = countNotifications + 1;
        end;
        
        % stable? J1 - J9 Jacobian matrix entries 

                J1 = -mu_par*(k_par/((k_par+N_LLM2_vector(j))^2))*P_LLM2_vector(j)-m_r_par/M(j);
                J2 = (-mu_par*u_LLM2_vector(j)+gamma_par*Phi_P_par)+(1-beta_par)*((2*g_par^2*eps_par*P_LLM2_vector(j))/((g_par+eps_par*(P_LLM2_vector(j))^2)^2))*Omega_par*Z_LLM2_vector(j);
                J3 = ((1-beta_par)*G_LLM2_vector(j)+Phi_Z_par)*Omega_par; 
                J4 = mu_par*(k_par/((k_par+N_LLM2_vector(j))^2))*P_LLM2_vector(j);
                J5 = mu_par*u_LLM2_vector(j)-Phi_P_par-m_r_par/M(j)-((2*g_par^2*eps*P_LLM2_vector(j))/((g_par+eps_par*(P_LLM2_vector(j))^2)^2))*Z_LLM2_vector(j);
                J6 = -G_LLM2_vector(j);
                J7 = 0;
                J8 = beta_par*((2*g_par^2*eps*P_LLM2_vector(j))/(((g_par+eps_par*(P_LLM2_vector(j))^2)^2)))*Z_LLM2_vector(j);
                J9 = beta_par*G_LLM2_vector(j)-Phi_Z_par-m_r_par/M(j);


                %Jacobian matrix
                J = [J1 J2 J3; J4 J5 J6; J7 J8 J9];

                % Eigen values of J
                lambda = eig(J);

                % test eigenvalues if equilibrium is stable. 2 is default
                % value for "not yet tested" 
                stable_LLM2_vector(j) = 2;

                for lam = 1:3
                    if abs(real(lambda(lam))) < eig_near_zero
                        stable_LLM2_vector(j) = 0;     % no conclusion possible
                    end;
                end;

                if stable_LLM2_vector(j) == 2
                   for lamb = 1:3 
                       if real(lambda(lamb)) > 0 
                           stable_LLM2_vector(j)= 1;      % unstable 
                       end;
                   end;
                end;

                if stable_LLM2_vector(j) == 2
                   stable_LLM2_vector(j) = -1;             % stable
                end;
        %%%%%%%%%%%%%%%%% put results into tables 
       if ND == ND_vector(1)

                    I_1_LLM2_table(3,j) = N_LLM2_vector(j);
                    I_1_LLM2_table(4,j) = dotN_LLM2_vector(j);
                    I_1_LLM2_table(5,j) = P_LLM2_vector(j);
                    I_1_LLM2_table(6,j) = dotP_LLM2_vector(j);
                    I_1_LLM2_table(7,j) = Z_LLM2_vector(j);
                    I_1_LLM2_table(8,j) = dotZ_LLM2_vector(j);                        
                    I_1_LLM2_table(9,j) = assumption1_LLM2_vector(j);
                    I_1_LLM2_table(10,j) = assumption2_LLM2_vector(j);
                    I_1_LLM2_table(11,j) = u_LLM2_vector(j);
                    I_1_LLM2_table(12,j) = G_LLM2_vector(j);
                    I_1_LLM2_table(13,j) = FN_LLM2_vector(j);
                    I_1_LLM2_table(14,j) = stable_LLM2_vector(j);

       else
           if ND == ND_vector(2)

                        I_8_LLM2_table(3,j) = N_LLM2_vector(j);
                        I_8_LLM2_table(4,j) = dotN_LLM2_vector(j);
                        I_8_LLM2_table(5,j) = P_LLM2_vector(j);
                        I_8_LLM2_table(6,j) = dotP_LLM2_vector(j);
                        I_8_LLM2_table(7,j) = Z_LLM2_vector(j);
                        I_8_LLM2_table(8,j) = dotZ_LLM2_vector(j);
                        I_8_LLM2_table(9,j) = assumption1_LLM2_vector(j);
                        I_8_LLM2_table(10,j) = assumption2_LLM2_vector(j);
                        I_8_LLM2_table(11,j) = u_LLM2_vector(j);
                        I_8_LLM2_table(12,j) = G_LLM2_vector(j);
                        I_8_LLM2_table(13,j) = FN_LLM2_vector(j);
                        I_8_LLM2_table(14,j) = stable_LLM2_vector(j);


           else         %  ND==3                       

                        I_16_LLM2_table(3,j) = N_LLM2_vector(j);
                        I_16_LLM2_table(4,j) = dotN_LLM2_vector(j);
                        I_16_LLM2_table(5,j) = P_LLM2_vector(j);
                        I_16_LLM2_table(6,j) = dotP_LLM2_vector(j);
                        I_16_LLM2_table(7,j) = Z_LLM2_vector(j);
                        I_16_LLM2_table(8,j) = dotZ_LLM2_vector(j);
                        I_16_LLM2_table(9,j) = assumption1_LLM2_vector(j);
                        I_16_LLM2_table(10,j) = assumption2_LLM2_vector(j);
                        I_16_LLM2_table(11,j) = u_LLM2_vector(j);
                        I_16_LLM2_table(12,j) = G_LLM2_vector(j);
                        I_16_LLM2_table(13,j) = FN_LLM2_vector(j);
                        I_16_LLM2_table(14,j) = stable_LLM2_vector(j);


           end;         %  if ND=2
       end;         %  if ND=1
    end;         % j = 1:400
end;        % loop ND
disp(['count Notifications ', num2str(countNotifications)]);

%%%%%%%%%%%%%%%%%%%%%
% create figures! 
%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(2,2,1)
plot(I_1_LLM2_table(1,:),I_1_LLM2_table(9,:),'k');
refline(0,0);
title('LLM2, Assumption1');
xlabel('M');
xlim([1,400]);
ylabel('Assumption1');


subplot(2,2,2)
hold on;
scatter(I_1_LLM2_table(1,:),I_1_LLM2_table(9,:),'k');
refline(0,0);
title('LLM2, Assumption 1, Detail');
xlabel('M');
xlim([1,5]);
ylabel('Assumption 1, Detail');
hold off;

subplot(2,2,3)
plot(I_1_LLM2_table(1,:),I_1_LLM2_table(10,:),'k');
refline(0,0);
title('LLM2, Assumption 2');
xlabel('M');
xlim([1,400]);
ylabel('Assumption 2');

figure(2)
subplot(3,2,1)
hold on;
plot(I_1_LLM2_table(1,:),I_1_LLM2_table(3,:),'b');
title('LLM2, N');
xlabel('M');
xlim([1,400]);
ylabel('N');

subplot(3,2,2)
hold on;
plot(I_1_LLM2_table(1,:),I_1_LLM2_table(3,:),'b');
title('LLM2, N (Detail)');
xlabel('M');
xlim([1,20]);
ylabel('N');
ylim([0,0.15]);

subplot(3,2,3)
hold on;
plot(I_1_LLM2_table(1,:),I_1_LLM2_table(5,:),'b');
plot(I_8_LLM2_table(1,:),I_8_LLM2_table(5,:),'r');
plot(I_16_LLM2_table(1,:),I_16_LLM2_table(5,:),'g');
title('LLM2, P');
xlabel('M');
xlim([1,400]);
ylabel('P');

subplot(3,2,5)
hold on;
plot(I_1_LLM2_table(3,:),I_1_LLM2_table(5,:),'b');
plot(I_8_LLM2_table(3,:),I_8_LLM2_table(5,:),'r');
plot(I_16_LLM2_table(3,:),I_16_LLM2_table(5,:),'g');
title('LLM2, N, P');
xlabel('N');
ylabel('P');
legend ('N_D = 1.5', 'N_D = 8', 'N_D = 16','Location','NorthEast');

subplot(3,2,6)
hold on;
plot(I_1_LLM2_table(3,:),I_1_LLM2_table(5,:),'b');
plot(I_8_LLM2_table(3,:),I_8_LLM2_table(5,:),'r');
plot(I_16_LLM2_table(3,:),I_16_LLM2_table(5,:),'g');
title('LLM2, N, P, Detail');
xlim([0,0.025]);
xlabel('N');
ylabel('P');

figure(3)

subplot(3,2,1)
scatter(I_1_LLM2_table(1,:),I_1_LLM2_table(14,:),'r');
title('LLM2, Exp.1, ND=1.5, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(3,2,2)
scatter(I_1_LLM2_table(1,:),I_1_LLM2_table(14,:),'r');
title('LLM2, Exp.1, ND=1.5, Stability, Detail');
xlim([0,5]);
xlabel('M');
xlim([1,5]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(3,2,3)
scatter(I_8_LLM2_table(1,:),I_8_LLM2_table(14,:),'r');
title('LLM2, Exp.1, ND=8, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(3,2,4)
scatter(I_8_LLM2_table(1,:),I_8_LLM2_table(14,:),'r');
title('LLM2, Exp.1, ND=8, Stability, Detail');
xlim([0,5]);
xlabel('M');
xlim([1,5]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(3,2,5)
scatter(I_16_LLM2_table(1,:),I_16_LLM2_table(14,:),'r');
title('LLM2, Exp.1, ND=16, Stability');
xlabel('M');
xlim([1,400]);
ylim([-1.2,1.2]);
ylabel('Stability');

subplot(3,2,6)
scatter(I_16_LLM2_table(1,:),I_16_LLM2_table(14,:),'r');
title('LLM2, Exp.1, ND=16, Stability, Detail');
xlim([0,5]);
xlabel('M');
xlim([1,5]);
ylim([-1.2,1.2]);
ylabel('Stability');

disp 'ready';
diary off;

