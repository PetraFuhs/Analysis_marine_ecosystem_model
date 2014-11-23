% QQM2, QQM3, QQM4

% Calculating Equilibria NPZ Model QQM for M = 1:400
% checks stability
% QQM ND = 1.5 or 8 or 16; 2 Parameter sets (Exp.II,III)
% Diploma Thesis P.Fuhs (CAU)
% Matlab Student Version V 7.10.0 (R2010a)

clear all;
clc;
clf;

delete('QQM_diary');
diary on; 
diary('QQM_diary');
countNotifications = 0;

disp 'Calculation for QQM 2,3,4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = NaN(1,400);         % Layer depth from 1 to 400 m, grid = 1m 

for j= 1:400
    M(j)=j;
end; 

% All units as in thesis.
% All domains as in thesis.

% Vector ND values: vector of three N_D values

ND_vector = [1.5, 8, 16];

% Parameter sets: Exp. II, III; two values of parameters per vector, experiment numbers as in [Schartau]

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

% how close to zero must an eigenvalue be to consider it zero? 
eig_near_zero = 10^-5;

%how close to zero must a derivative be to consider it zero? 
deriv_near_zero = 0.0001;


% Create result tables for QQM2, QQM3, QQM4:
% 
% 6*3 = 18  tables: ND = 1.5,8,16; Exp. II,III => 6 cases
% in all 6 cases up to 3 solutions QQM2, QQM3, QQM4
%
% table names: ExpNumber_NDValue_QQMNumber_table, e.g.II_8_QQM3_table
% ND = 1 in table names refers to 1.5 actual value of ND

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
posa2 = 10;         % line 10: a2
posb2 = 11;         % line 11: b2
posc2 = 12;         % line 12: c2
postheta = 13;      % line 13: theta
posnumberreal = 14; % line 14: number of real roots 
posnumberusef = 15; % line 15: number of useful roots (real and nonnegative)
posu = 16;          % line 16: u
posG = 17;          % line 17: G
posFN = 18;         % line 18: FN
posstable = 19;     % line 19; stable -1=yes, 0=not determinded, 1=no, default: 2


II_1_QQM2_table = NaN(19,400);    II_1_QQM3_table = NaN(19,400);    II_1_QQM4_table = NaN(19,400);
II_8_QQM2_table = NaN(19,400);    II_8_QQM3_table = NaN(19,400);    II_8_QQM4_table = NaN(19,400);
II_16_QQM2_table = NaN(19,400);   II_16_QQM3_table = NaN(19,400);   II_16_QQM4_table = NaN(19,400);
III_1_QQM2_table = NaN(19,400);   III_1_QQM3_table = NaN(19,400);   III_1_QQM4_table = NaN(19,400);
III_8_QQM2_table = NaN(19,400);   III_8_QQM3_table = NaN(19,400);   III_8_QQM4_table = NaN(19,400);
III_16_QQM2_table = NaN(19,400);  III_16_QQM3_table = NaN(19,400);  III_16_QQM4_table = NaN(19,400);

% create auxiliary vector for results
aux_vector = NaN(19,1);   
     
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
    
    
    for exp = 2:3         % loop experiments

        % Display Model_StationaryPoint, ND, display experiment number
        disp ([' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m']);

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

        for j= 1:400         % loop layer depth M from 1m to 400 m, calculate terms and variables
            
            % fill in lines for M, ND, exp
            aux_vector(posM,1) = j;
            aux_vector(posND,1) = ND;
            aux_vector(posExp,1) = exp;
            % calculate auxiliary terms 
                
            % a2 
            aux_vector(posa2,1) =  (-M(j)/m_r_par)*((1/Phi_P_opt_par)*(mu_par*(-mu_par+Phi_P_par+m_r_par/M(j)+gamma_par*Phi_P_par)+gamma_par*Phi_P_par*(-Phi_P_par-m_r_par/M(j)))+m_r_par/M(j)*(ND-2*k_par));

            % b2
            aux_vector(posb2,1) = (-M(j)/m_r_par)*((k_par/Phi_P_opt_par)*(Phi_P_par*(mu_par+gamma_par*mu_par-2*gamma_par*Phi_P_par-2*gamma_par*(m_r_par/M(j)))+mu_par*m_r_par/M(j))+(k_par*m_r_par/M(j))*(2*ND-k_par));
            
            % c2 
            aux_vector(posc2,1) = ((gamma_par*Phi_P_par*((k_par)^2)*M(j))/(Phi_P_opt_par*m_r_par))*(Phi_P_par+m_r_par/M(j))-ND*((k_par)^2);

            % theta 
            aux_vector(postheta,1) = 18*aux_vector(posa2,1)*aux_vector(posb2,1)*aux_vector(posc2,1) + (aux_vector(posa2,1))^2*(aux_vector(posb2,1))^2 - 27*(aux_vector(posc2,1))^2 - 4*(aux_vector(posb2,1))^3 - 4*(aux_vector(posa2,1))^3*aux_vector(posc2,1);

            % find out number of real roots (aux_vector(posnumberreal,1))
            %("Diskriminantensatz", [Dörrie], theta = line 1)
            % "useful roots" in this case are real and within domain 
            % => pointer 2 might be smaller than
            % aux_vector(posnumberreal,1)
            
            if aux_vector(postheta,1) > 0 
                aux_vector(posnumberreal,1) = 3;
            else         % not 3 real roots, theta not greater zero
                if aux_vector(postheta,1) < 0
                    aux_vector(posnumberreal,1) = 1;
                else         % theta = 0 : one root is double!
                    aux_vector(posnumberreal,1) = 2;
                end         % not 3 real roots
            end;         % end theta -> aux_vector(posnumberreal,1)

            % calculate variables 

            % have MATLAB find roots of N numerically; p = polynom N

            % p = 1 N^3 + a2 N^2 + b2 N + c2
            p = [1 aux_vector(posa2,1) aux_vector(posb2,1) aux_vector(posc2,1)];
            
            % MATLAB finds roots via eigen values 
            rootsN_QQM_vector = roots(p);

            % sort roots; real roots , nonnegative, are useful; rest not usable; usable roots first in vector!
            sort_rootsN_QQM_vector = rootsN_QQM_vector;
            pointer1 = 1;         % pointer1: next position for useful root (start with position 1)
            pointer2 = 3;         % pointer2: next position for not useful root (start with position 3)

            condition1 = NaN;
            condition2 = NaN;
            
            for w = 1:3         % sort the roots of N into a sorted vector: useful first 
                
                % conditions for acceptance of root
                if abs(imag(rootsN_QQM_vector(w))) < eig_near_zero         % looking for real root
                    condition1 = 1;
                else         % root not real
                    condition1 = 0; 
                end;         % looking for real root

                if  real(rootsN_QQM_vector(w))>= 0         % looking for nonnegative root
                    condition2 = 1;
                else         % root negative 
                    condition2 = 0; 
                end;         % looking for nonnegative root 
                
                condition = condition1*condition2; 
               
                if  condition == 1         % real and nonnegative
                    
                    sort_rootsN_QQM_vector(pointer1)=rootsN_QQM_vector(w);
                    pointer1 = pointer1+1;
                   
                else         % not real, or negative 
                    sort_rootsN_QQM_vector(pointer2)=rootsN_QQM_vector(w);
                    pointer2 = pointer2-1;
                    
                end;         % root real and nonnegative 
            end;         % w = 1:3, sorting roots, useful first, not useful last 
            
            % real values of N (within domain) -> find P and construct
            % (N_QQM, P_QQM, 0), u, FN, G, 
            % and dotN, dotP, dotZ (if proper layer depth M then dotVar almost zero)
            % and fill in table / tables 
            
            % number of usable roots in pointer2
              
            if pointer2 == 0         % pointer2 = 0 
               disp ([' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  all roots are complex or out of domain']);
               countNotifications =  countNotifications + 1;
               return;
               
            else         % pointer2 > 0 % pointer2 not equal 0
               aux_vector(posnumberusef,1) = pointer2;
            end;         % pointer2 = 0?
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % start calculation N,P,Z, dotN, dotP, dotZ, u, G, FN
            % for QQM2
            % if pointer2 >1, same calculation for QQM3 
            % if pointer2 >2, same calculation for QQM4 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % N : aux_vector(posN,1)
            
            for pointnumber = 1:aux_vector(posnumberusef,1)   %loop pointnumber 1:# of useful roots
                
                aux_vector(posN,1) = sort_rootsN_QQM_vector(pointnumber);
                numberofequil = pointnumber + 1;
                % N > ND
                if aux_vector(posN,1) > ND
                    disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m , |  N is greater than ND']);
                    countNotifications =  countNotifications + 1;
                end;         % N > ND

                % u : aux_vector(posu,1) 
                aux_vector(posu,1) = aux_vector(posN,1)/(aux_vector(posN,1)+k_par);

                % u >= 1 
                if aux_vector(posu,1) >= 1
                   disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u is greater or equal 1']);
                   countNotifications =  countNotifications + 1;
                end;         
                % u >= 1 
                if aux_vector(posu,1) < 0 
                   disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  u is negative']);
                   countNotifications = countNotifications + 1;
                end;

                % FN : aux_vector(posFN,1)
                aux_vector(posFN,1) = (m_r_par/M(j))*(ND-aux_vector(posN,1));

                % FN < 0
                if aux_vector(posFN,1) < 0 
                    disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  F_N is negative']);
                    countNotifications =  countNotifications + 1;
                end;         % FN negative

                % FN too big? 
                if aux_vector(posFN,1) > 60
                    disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  F_N > 60']);
                    countNotifications = countNotifications + 1;
                end;        % FN too big 

                % P : aux_vector(posP,1)
                aux_vector(posP,1) = (1/Phi_P_opt_par)*(mu_par*aux_vector(posu,1)-Phi_P_par-(m_r_par/M(j)));

                % P < 0
                if aux_vector(posP,1) < 0 
                    disp(['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  P is negative']);
                    countNotifications = countNotifications + 1;
                end;        % P negative

                % G: aux_vector(posG,1)
                aux_vector(posG,1) = (g_par*eps_par*(aux_vector(posP,1))^2)/(g_par+eps_par*(aux_vector(posP,1))^2);

                % G negative (e.g. parameter from wrong domain) 
                if aux_vector(posG,1) < 0
                   disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | Notification. G is negative']);
                   countNotifications =  countNotifications + 1;
                end;         % G negative

                % G greater g
                if aux_vector(posG,1) >= g_par
                   disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  G is greater than parameter g']);
                   countNotifications =  countNotifications + 1;
                end;         % G greater g

                % Z: aux_vector(posZ,1)
                aux_vector(posZ,1)= 0;

                % dotN: aux_vector(posdotN,1)
                aux_vector(posdotN,1)= (-mu_par*aux_vector(posu,1)+gamma_par*Phi_P_par)*aux_vector(posP,1)+((1-beta_par)*aux_vector(posG,1)+Phi_Z_par)*Omega_par*aux_vector(posZ,1)+aux_vector(posFN,1);

                % dotN not almost zero
                if abs(aux_vector(posdotN,1)) > deriv_near_zero
                    disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, | N might not be stationary at this depth']);
                    countNotifications = countNotifications +1;
                end;         % dotN not zero

                % dotP: aux_vector(posdotP,1)
                aux_vector(posdotP,1)= (mu_par*aux_vector(posu,1)-Phi_P_par-(m_r_par/M(j)))*aux_vector(posP,1)-aux_vector(posG,1)*aux_vector(posZ,1)-Phi_P_opt_par*(aux_vector(posP,1))^2;

                % dotP not almost zero
                if abs(aux_vector(posdotP,1)) > deriv_near_zero
                    disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  P might not be stationary at this depth']);
                    countNotifications = countNotifications +1;
                end;

                % dotZ = aux_vector(posdotZ,1) automatically zero here but
                % formula given for completeness 
                aux_vector(posdotZ,1)= (beta_par*aux_vector(posG,1)-Phi_Z_par-(m_r_par/M(j)))*aux_vector(posZ,1)-Phi_Z_opt_par*(aux_vector(posZ,1))^2;

                % dotZ not almost zero
                if abs(aux_vector(posdotZ,1)) > deriv_near_zero
                    disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Z might not be stationary at this depth']);
                    countNotifications = countNotifications +1;
                end;

                % test for positivity in case beta > 1
                if ((1-beta_par)*aux_vector(posG,1) + Phi_Z_par)*Omega_par*aux_vector(posZ,1) <0
                   if  gamma_par*Phi_P_par*aux_vector(posP,1) + m_r_par/M(j)*ND < abs(((1-beta_par)*aux_vector(posG,1) + Phi_Z_par)*Omega_par*aux_vector(posZ,1))
                    disp (['QQM',num2str(numberofequil), ' ND: ',num2str(ND),' exp: ',num2str(exp),' depth: ',num2str(j),' m, |  Positivity is violated']);
                   countNotifications = countNotifications + 1;
                   end;
                end;

                % stable? J1 - J9 Jacobian matrix entries 

                J1 = -mu_par*(k_par/((k_par+aux_vector(posN,1))^2))*aux_vector(posP,1)-m_r_par/M(j);
                J2 = (-mu_par*aux_vector(posu,1)+gamma_par*Phi_P_par)+(1-beta_par)*(2*g_par^2*eps_par*aux_vector(posP,1)/((g_par+eps_par*aux_vector(posP,1)^2)^2))*Omega_par*aux_vector(posZ,1);
                J3 = ((1-beta_par)*aux_vector(posG,1)+Phi_Z_par)*Omega_par; 
                J4 = mu_par*(k_par/((k_par+aux_vector(posN,1))^2))*aux_vector(posP,1);
                J5 = mu_par*aux_vector(posu,1)-Phi_P_par-m_r_par/M(j)-(2*g_par^2*eps*aux_vector(posP,1)/((g_par+eps_par*aux_vector(posP,1)^2)^2))*aux_vector(posZ,1)-2*Phi_P_opt_par*aux_vector(posP,1);
                J6 = -aux_vector(posG,1);
                J7 = 0;
                J8 = beta_par*(2*g_par^2*eps*aux_vector(posP,1)/((g_par+eps_par*aux_vector(posP,1)^2)^2))*aux_vector(posZ,1);
                J9 = beta_par*aux_vector(posG,1)-Phi_Z_par-m_r_par/M(j)-2*Phi_Z_opt_par*aux_vector(posZ,1);

                %Jacobian matrix
                J = [J1 J2 J3; J4 J5 J6; J7 J8 J9];

                % Eigen values of J
                lambda = eig(J);

                % test eigenvalues if equilibrium is stable. 2 is default
                % value for "not yet tested" 
                aux_vector(posstable,1) = 2;

                for lam = 1:3
                    if abs(real(lambda(lam))) < eig_near_zero
                        aux_vector(posstable,1) = 0;     % no conclusion possible
                    end;
                end;

                if aux_vector(posstable,1) == 2
                   for lamb = 1:3 
                       if real(lambda(lamb)) > 0 
                           aux_vector(posstable,1)= 1;      % unstable 
                       end;
                   end;
                end;

                if aux_vector(posstable,1) == 2
                   aux_vector(posstable,1) = -1;             % stable
                end;
           
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End QQM 2,3,4 calculating
           
                %%%%%%%%%%%% putting results into table: 
                if pointnumber== 1 
                    if ND == ND_vector(1)

                        if exp == 2
                            II_1_QQM2_table(:,j) = aux_vector;
                        else         % exp = 3
                            III_1_QQM2_table(:,j) = aux_vector;
                        end;         % exp 

                    else         % ND not 1.5

                        if ND == ND_vector(2)
                            if exp == 2
                                II_8_QQM2_table(:,j) = aux_vector;
                            else         % exp = 3
                                III_8_QQM2_table(:,j) = aux_vector;
                            end;         % exp 

                        else         % ND == ND_vector(3)
                            if exp == 2
                                II_16_QQM2_table(:,j) = aux_vector;
                            else         % exp = 3
                                III_16_QQM2_table(:,j) = aux_vector;
                            end;         % exp 
                        end;         % ND not 1.5
                    end;         % ND = 1.5
                else 
                    if pointnumber== 2 
                     
                       if ND == ND_vector(1)

                           if exp == 2
                               II_1_QQM3_table(:,j) = aux_vector;
                           else         % exp = 3
                               III_1_QQM3_table(:,j) = aux_vector;
                           end;         % exp 

                       else         % ND not 1.5

                           if ND == ND_vector(2)
                               if exp == 2
                                   II_8_QQM3_table(:,j) = aux_vector;
                               else         % exp = 3
                                   III_8_QQM3_table(:,j) = aux_vector;
                               end;         % exp 

                           else         % ND == ND_vector(3)
                               if exp == 2
                                   II_16_QQM3_table(:,j) = aux_vector;
                               else         % exp = 3
                                   III_16_QQM3_table(:,j) = aux_vector;
                               end;         % exp 
                           end;         % ND not 1.5
                       end;         % ND = 1.5
                    else   

                       if ND == ND_vector(1)

                           if exp == 2
                               II_1_QQM4_table(:,j) = aux_vector;
                           else         % exp = 3
                               III_1_QQM4_table(:,j) = aux_vector;
                           end;         % exp 

                       else         % ND not 1.5

                           if ND == ND_vector(2)
                               if exp == 2
                                   II_8_QQM4_table(:,j) = aux_vector;
                               else         % exp = 3
                                   III_8_QQM4_table(:,j) = aux_vector;
                               end;         % exp 

                           else         % ND == ND_vector(3)
                               if exp == 2
                                   II_16_QQM4_table(:,j) = aux_vector;
                               else         % exp = 3
                                   III_16_QQM4_table(:,j) = aux_vector;
                               end;         % exp 
                           end;         % ND not 1.5
                        end;         % ND = 1.5
                    end;    % pointnumber= 2  
                end;     % pointnumber= 1
            end;         % loop pointnumber = 1:pointer2
        end;         % loop layer depth
    end;        % loop experiments
end;        % loop ND


               

disp ' Notifications : '
disp(countNotifications);

%%%%%%%%%%%%%%%%
% QQM2, Exp. II
% Plots: first plots are N to M, P to M, each with 3 values of ND
% second: N,P to M, with each one value of ND (3 times) 
% third: P to N with 3 ND 

figure(1)

subplot (2,2,1) 
hold on; 
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(4,:),'b');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(4,:),'r');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(4,:),'g');
title('QQM2, Exp.II, N');
xlabel('M');
xlim([1,400]);
ylabel('N');
legend('N_D = 1.5','N_D = 8','N_D = 16');

subplot (2,2,2) 
hold on; 
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(4,:),'b');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(4,:),'r');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(4,:),'g');
title('QQM2, Exp.II, N, Detail');
xlabel('M');
xlim([1,20]);
ylabel('N');

subplot (2,2,3) 
hold on; 
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(6,:),'b');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(6,:),'r');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(6,:),'g');
title('QQM2, Exp.II, P');
xlabel('M');
xlim([1,400]);
ylabel('P');
ylim([-1,20]);

subplot (2,2,4) 
hold on; 
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(6,:),'b');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(6,:),'r');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(6,:),'g');
title('QQM2, Exp.II, P, Detail');
xlabel('M');
xlim([1,30]);
ylabel('P');
ylim([-1,20]);


figure(2)
 
subplot(1,2,1)
hold on;
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(4,:),':b');
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(6,:),'--b');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(4,:),':r');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(6,:),'--r');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(4,:),':g');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(6,:),'--g');
title('QQM2, Exp.II, N and P on M');
xlabel('M');
xlim([1,400]);
ylabel('N,P');
legend('N_D = 1.5, N','N_D = 1.5, P','N_D = 8, N','N_D = 8, P','N_D = 16, N','N_D = 16, P');

subplot(1,2,2)
hold on;
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(4,:),':b');
plot(II_1_QQM2_table(1,:),II_1_QQM2_table(6,:),'--b');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(4,:),':r');
plot(II_8_QQM2_table(1,:),II_8_QQM2_table(6,:),'--r');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(4,:),':g');
plot(II_16_QQM2_table(1,:),II_16_QQM2_table(6,:),'--g');
title('QQM2, Exp.II, N and P on M');
xlabel('M');
xlim([1,10]);
ylabel('N,P');



figure(3)

hold on;
plot(II_1_QQM2_table(4,:),II_1_QQM2_table(6,:),'b');
plot(II_8_QQM2_table(4,:),II_8_QQM2_table(6,:),'r');
plot(II_16_QQM2_table(4,:),II_16_QQM2_table(6,:),'g');
title('QQM2, Exp.II, P on N');
xlabel('N');
ylabel('P');
legend('N_D = 1.5','N_D = 8','N_D = 16');

%%%%%%%%%%%%%%%%
% QQM2, Exp. III
% Plots: first plots are N to M, P to M, each with 3 values of ND
% second: N,P to M, with each one value of ND (3 times) 
% third: P to N with 3 ND 

figure(4)

subplot (2,2,1) 
hold on; 
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(4,:),'b');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(4,:),'r');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(4,:),'g');
title('QQM2, Exp.III, N');
xlabel('M');
xlim([1,400]);
ylabel('N');
legend('N_D = 1.5','N_D = 8','N_D = 16');

subplot (2,2,2) 
hold on; 
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(4,:),'b');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(4,:),'r');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(4,:),'g');
title('QQM2, Exp.III, N, Detail');
xlabel('M');
xlim([1,20]);
ylabel('N');

subplot (2,2,3) 
hold on; 
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(6,:),'b');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(6,:),'r');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(6,:),'g');
title('QQM2, Exp.III, P');
xlabel('M');
xlim([1,400]);
ylabel('P');

subplot (2,2,4) 
hold on; 
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(6,:),'b');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(6,:),'r');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(6,:),'g');
title('QQM2, Exp.III, P, Detail');
xlabel('M');
xlim([1,20]);
ylabel('P');

figure(5)
 
subplot(1,2,1)
hold on;
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(4,:),':b');
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(6,:),'--b');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(4,:),':r');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(6,:),'--r');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(4,:),':g');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(6,:),'--g');
title('QQM2, Exp.III, N and P on M');
xlabel('M');
xlim([1,400]);
ylabel('N,P');
legend('N_D = 1.5, N','N_D = 1.5, P','N_D = 8, N','N_D = 8, P','N_D = 16, N','N_D = 16, P');

subplot(1,2,2)
hold on;
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(4,:),':b');
plot(III_1_QQM2_table(1,:),III_1_QQM2_table(6,:),'--b');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(4,:),':r');
plot(III_8_QQM2_table(1,:),III_8_QQM2_table(6,:),'--r');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(4,:),':g');
plot(III_16_QQM2_table(1,:),III_16_QQM2_table(6,:),'--g');
title('QQM2, Exp.III, N and P on M');
xlabel('M');
xlim([1,10]);
ylabel('N,P');



figure(6)

hold on;
plot(III_1_QQM2_table(4,:),III_1_QQM2_table(6,:),'b');
plot(III_8_QQM2_table(4,:),III_8_QQM2_table(6,:),'r');
plot(III_16_QQM2_table(4,:),III_16_QQM2_table(6,:),'g');
title('QQM2, Exp.III, P on N');
xlabel('N');
ylabel('P');
legend('N_D = 1.5','N_D = 8','N_D = 16');

figure(7) %(stability exp.2) 

subplot(2,2,1)
hold on;
scatter(II_1_QQM2_table(posM,:),II_1_QQM2_table(posstable,:),'r');
title('QQM2, Exp.II, ND=1.5, Stability');
xlabel('M');
xlim([1,400]);
ylabel('Stability');
ylim([-1.2,1.2]);
legend('-1=stable;+1=unstable;0=not determined');

subplot(2,2,2)
hold on;
scatter(II_1_QQM2_table(posM,:),II_1_QQM2_table(posstable,:),'r');
title('QQM2, Exp.II, ND=1.5, Stability, Detail');
xlabel('M');
xlim([290,300]);
ylabel('Stability');
ylim([-1.2,1.2]);

subplot(2,2,3)
hold on;
scatter(II_8_QQM2_table(posM,:),II_8_QQM2_table(posstable,:),'r');
title('QQM2, Exp.II, ND=8, Stability');
xlabel('M');
xlim([1,400]);
ylabel('Stability');
ylim([-1.2,1.2]);

subplot(2,2,4)
hold on;
scatter(II_16_QQM2_table(posM,:),II_16_QQM2_table(posstable,:),'r');
title('QQM2, Exp.II, ND=16, Stability');
xlabel('M');
xlim([1,400]);
ylabel('Stability');
ylim([-1.2,1.2]);

figure(8) % stability exp.III
subplot(2,2,1)
hold on;
scatter(III_1_QQM2_table(posM,:),III_1_QQM2_table(posstable,:),'r');
title('QQM2, Exp.III, ND=1.5, Stability');
xlabel('M');
xlim([1,400]);
ylabel('Stability');
ylim([-1.2,1.2]);
 
subplot(2,2,2)
hold on;
scatter(III_1_QQM2_table(posM,:),III_1_QQM2_table(posstable,:),'r');
title('QQM2, Exp.III, ND=1.5, Stability, Detail');
xlabel('M');
xlim([1,5]);
ylabel('Stability');
ylim([-1.2,1.2]);

subplot(2,2,3)
hold on;
scatter(III_8_QQM2_table(posM,:),III_8_QQM2_table(posstable,:),'r');
title('QQM2, Exp.III, ND=8, Stability');
xlabel('M');
xlim([1,400]);
ylabel('Stability');
ylim([-1.2,1.2]);

subplot(2,2,4)
hold on;
scatter(III_16_QQM2_table(posM,:),III_16_QQM2_table(posstable,:),'r');
title('QQM2, Exp.III, ND=16, Stability');
xlabel('M');
xlim([1,400]);
ylabel('Stability');
ylim([-1.2,1.2]);


% make a nice plot of u, G, and FN:

% 
% subplot(1,3,1)
% plot(II_8_QQM2_table(4,:), II_8_QQM2_table(16,:));
% title('Exp.II, ND=8, QQM2');
% xlabel('N in mmolN/m^3');
% xlim([0.0203,0.2430]);
% ylabel('u, no unit');
% ylim([0,0.3]);
% 
% subplot(1,3,2)
% plot(II_8_QQM2_table(4,:), II_8_QQM2_table(18,:));
% title('Exp.II, ND=8, QQM2');
% xlabel('N in mmolN/m^3');
% xlim([0.0203,0.2430]);
% ylabel('F_N in mmolN/m^3');
% ylim([0,0.9]);
% 
% 
% subplot(1,3,3)
% plot(II_8_QQM2_table(6,:), II_8_QQM2_table(17,:));
% title('Exp.II, ND=8, QQM2');
% xlabel('P in mmolN/m^3');
% xlim([0.1100,1.9270]);
% ylabel('G in 1/d');
% ylim([0,0.95]);
% 



disp 'ready';
disp '********************************';
diary off

               
            
          
   