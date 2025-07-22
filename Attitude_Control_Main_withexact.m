format long
clear all; close all; clc;
%% Code to solve the attitude control problem
% Written by Bethany Hintz - 2 July 2025

% Initial Values
epsilon = 10;
N = 15;% Number of nodes
M = 3; % Number of elements
tf = 3.545; % Final time

% Initial Parameteres
Neq = 14;

phi = pi;

% Initial conditions for x_bar
% q10 q20 q30 q40 w10 w20 w30
x0 = [0, 0, 0, 1, 0, 0, 0];
%Final conditions for x_bar
% q1f q2f q3f q4f w1f w2f w3f
xf = [0,0,sin(phi/2),cos(phi/2),0,0,0];

BCtype = 'fixed';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'
BC = [x0,xf]';

initial_guess = [ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1);ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1);ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1);ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1)]; 

%Time guess
t1 = 1;
t2 = 3;
%% Approximation matrix
D = zeros(N,N,M);
phi = zeros(N,N,M);
phid = zeros(N,N,M);

for k = 1
    t0e            = (k-1)*tf/M; % Initial time of each segment (element)
    tfe            = t1;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end

for k = 2
    t0e            = t1; % Initial time of each segment (element)
    tfe            = t2;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end

for k = M
    t0e            = t2; % Initial time of each segment (element)
    tfe            = k*tf/M;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    %     tinterp(:,:,k) = linspace(t0e,tfe,5*n);
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
        %         phi_interp(i,:,k) = rbf0(c,te(:,i,k),tinterp(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end
 
%% CRBF Solution using fsolve
ip=0;
Jlocal_pattern = ones(Neq*N);
J_pattern = kron(eye(M),Jlocal_pattern);

J_local = zeros(Neq*N,Neq*N,M);
J    = zeros(Neq*N*M, Neq*N*M); % Local matrices are to be concatenated sequentially

ns = 0;
Elb   = zeros(1,M);
Erb   = Elb;
SLb   = zeros(Neq,1,M);

SRb   = SLb;
R = zeros(N,M);
sizeJ_local = size(J_local);
LJ_local    = sizeJ_local(1);

f = @(X) attitude_local(X,BC,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
options = optimoptions(@fsolve,'Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',false,'JacobPattern',J_pattern,'StepTolerance',1e-20,'FunctionTolerance',1e-20,'UseParallel',true,'FiniteDifferenceType','central');
[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);

x = reshape(xx,[Neq*N,M]);
Q1 = x(1:N,:);
Q2 = x(N+1:2*N,:);
Q3 = x(2*N+1:3*N,:);
Q4 = x(3*N+1:4*N,:);
W1 = x(4*N+1:5*N,:);
W2 = x(5*N+1:6*N,:);
W3 = x(6*N+1:7*N,:);
L1 = x(7*N+1:8*N,:);
L2 = x(8*N+1:9*N,:);
L3 = x(9*N+1:10*N,:);
L4 = x(10*N+1:11*N,:);
L5 = x(11*N+1:12*N,:);
L6 = x(12*N+1:13*N,:);
L7 = x(13*N+1:14*N,:);


ICs = [x0(1) x0(2) x0(3) x0(4) x0(5) x0(6) x0(7) L1(1) L2(1) L3(1) L4(1) L5(1) L6(1) L7(1)];
opts = odeset('RelTol',1e-14,'AbsTol',1e-20);
T1 = te(:);
q11 =0;
q22 =0;
q33 =0;
q44 =0;
q55 =0;
q66 =0;

[TTT, YYY] = ode45 (@(T1,x)attdynamics(x), T1, ICs);


%% Clean Data %%
T = te(:);

% [~, ind] = unique(T);
% duplicate_ind = flip(setdiff(1:size(T), ind));
% duplicate_value = T(duplicate_ind);
% for i = 1 : length(duplicate_ind)
%     T(duplicate_ind(i))  = [];
%     X1(duplicate_ind(i)) = [];
%     X2(duplicate_ind(i)) = [];
%     Y1(duplicate_ind(i)) = [];
%     Y2(duplicate_ind(i)) = [];
%     Z1(duplicate_ind(i)) = [];
%     Z2(duplicate_ind(i)) = [];
%     L1(duplicate_ind(i)) = [];
%     L2(duplicate_ind(i)) = [];
% 
% end
q1 = Q1(:);
q2 = Q2(:);
q3 = Q3(:);
q4 = Q4(:);
w1 = W1(:);
w2 = W2(:);
w3 = W3(:);

qq1 = YYY(:,1);
qq2 = YYY(:,2);
qq3 = YYY(:,3);
qq4 = YYY(:,4);
ww1 = YYY(:,5);
ww2 = YYY(:,6);
ww3 = YYY(:,7);
LL1 = YYY(:,8);
LL2 = YYY(:,9);
LL3 = YYY(:,10);
LL4 = YYY(:,11);
LL5 = YYY(:,12);
LL6 = YYY(:,13);
LL7 = YYY(:,14);

error_q1 = (q1 - qq1)./qq1;
error_q2 = (q2 - qq2)./qq2;
error_q3 = (q1 - qq3)./qq3;
error_q4 = (q2 - qq4)./qq4;
error_w1 = (w1 - ww1)./ww1;
error_w2 = (w2 - ww2)./ww2;
error_w3 = (w3 - ww3)./ww3;


abserror_q1 = (q1 - qq1);
abserror_q2 = (q2 - qq2);
abserror_q3 = (q3 - qq3);
abserror_q4 = (q4 - qq4);
abserror_w1 = (w1 - ww1);
abserror_w2 = (w2 - ww2);
abserror_w3 = (w3 - ww3);


%% Plot Results %%
figure(1)
plot(T,q1,'*','LineWidth',1.1)
hold on
plot(T,q2,'^','LineWidth',1.1)
plot(T,q3,'*','LineWidth',1.5)
plot(T,q4,'*','LineWidth',1.5)
plot(T,qq1,'-','LineWidth',1.1)
plot(T,qq2,'-','LineWidth',1.1)
plot(T,qq3,'-','LineWidth',1.5)
plot(T,qq4,'-','LineWidth',1.5)
xlabel('Time')
ylabel('Quarternion components')
title('Attitude Control - Quarternions')
legend('q1 (CRBF)','q2 (CRBF)','q3 (CRBF)','q4 (CRBF)','q1 (Exact)', 'q2 (Exact)', 'q3 (Exact)', 'q4 (Exact)')

figure(2)
plot(T,w1,'*','LineWidth',1.1)
hold on
plot(T,w2,'^','LineWidth',1.1)
plot(T,w3,'*','LineWidth',1.5)
plot(T,ww1,'-','LineWidth',1.1)
plot(T,ww2,'-','LineWidth',1.1)
plot(T,ww3,'-','LineWidth',1.5)
xlabel('Time')
ylabel('Angular Velocity')
title('Attitude Control - Angular Velocities')
legend('\omega_1 (CRBF)','\omega_2 (CRBF)','\omega_3 (CRBF)','\omega_1 (Exact)', '\omega_2 (Exact)', '\omega_3 (Exact)')

%Orientation angles
alpha = acos(1-2*q2.^2-2*q3.^2);
beta = acos(1-2*q3.^2-2*q1.^2);
gamma = acos(1-2*q1.^2-2*q2.^2);
alpha = rad2deg(alpha);
beta = rad2deg(beta);
gamma = rad2deg(gamma);

alphaa = acos(1-2*qq2.^2-2*qq3.^2);
betaa = acos(1-2*qq3.^2-2*qq1.^2);
gammaa = acos(1-2*qq1.^2-2*qq2.^2);
alphaa = rad2deg(alphaa);
betaa = rad2deg(betaa);
gammaa = rad2deg(gammaa);

abserror_alpha = (alpha - alphaa);
abserror_beta = (beta - betaa);
abserror_gamma = (gamma - gammaa);

error_alpha = (alpha - alphaa)./alphaa;
error_beta = (beta - betaa)./betaa;
error_gamma = (gamma - gammaa)./gammaa;

figure(3)
plot(T,alpha,'^','LineWidth',1.5)
hold on
plot(T,beta,'*','LineWidth',0.5)
plot(T,gamma,'*','LineWidth',1.5)
plot(T,alphaa,'-','LineWidth',1.5)
plot(T,betaa,'-','LineWidth',0.5)
plot(T,gammaa,'-','LineWidth',1.5)
xlabel('Time')
ylabel('Orientation angle (deg)')
title('Attitude Control - Orientation Angles')
legend('\alpha (CRBF)','\beta (CRBF)','\gamma (CRBF)', '\alpha (Exact)','\beta (Exact)','\gamma (Exact)')

figure(4)
plot(T,L5(:),'*','LineWidth',0.5)
hold on
plot(T,L6(:),'^','LineWidth',1.5)
plot(T,L7(:),'*','LineWidth',1.5)
plot(T,LL5(:),'-','LineWidth',0.5)
plot(T,LL6(:),'-','LineWidth',1.5)
plot(T,LL7(:),'-','LineWidth',1.5)
xlabel('Time')
ylabel('Control Components along Body Axis')
title('Attitude Control - Control Components')
legend('u1 (CRBF)','u2 (CRBF)','u3 (CRBF)', 'u1 (Exact)','u2 (Exact)','u3 (Exact)')


% figure;
% plot(T(2:end),abs(error_q1(2:end)))
% hold on
% plot(T(2:end),abs(error_q2(2:end)))
% plot(T(2:end),abs(error_q3(2:end)))
% plot(T(2:end),abs(error_q4(2:end)))
% xlabel('Time')
% ylabel('Relative Error')
% title('Quarternion Relative Error')
% legend('q1', 'q2', 'q3', 'q4')
% 
% figure;
% plot(T(2:end),abs(error_w1(2:end)))
% hold on
% plot(T(2:end),abs(error_w2(2:end)))
% plot(T(2:end),abs(error_w3(2:end)))
% xlabel('Time')
% ylabel('Relative Error')
% title('Angular Velocity Component Relative Error')
% legend('\omega_1', '\omega_2', '\omega_3')

figure;
plot(T,abs(abserror_q1),'-*','LineWidth',1.1)
hold on
plot(T,abs(abserror_q2),'-^','LineWidth',1.1)
plot(T,abs(abserror_q3),'-*','LineWidth',1.1)
plot(T,abs(abserror_q4),'-*','LineWidth',1.1)
xlabel('Time')
ylabel('Absolute Error')
title('Quarternion Absolute Error')
legend('q1', 'q2', 'q3', 'q4')

figure;
plot(T,abs(abserror_w1),'-*','LineWidth',1.1)
hold on
plot(T,abs(abserror_w2),'-^','LineWidth',1.1)
plot(T,abs(abserror_w3),'-*','LineWidth',1.1)
xlabel('Time')
ylabel('Absolute Error')
title('Angular Velocity Component Absolute Error')
legend('\omega_1', '\omega_2', '\omega_3')

figure;
plot(T,abs(abserror_alpha),'-*','LineWidth',1.1)
hold on
plot(T,abs(abserror_beta),'-^','LineWidth',1.1)
plot(T,abs(abserror_gamma),'-*','LineWidth',1.1)
xlabel('Time')
ylabel('Absolute Error (deg)')
title('Orientation Angles Absolute Error')
legend('\alpha', '\beta', '\gamma')

% figure;
% plot(T(2:end),abs(error_alpha(2:end)))
% hold on
% plot(T(2:end),abs(error_beta(2:end)))
% plot(T(2:end),abs(error_gamma(2:end)))
% xlabel('Time')
% ylabel('Relative Error')
% title('Angular Velocity Component Relative Error')
% legend('\alpha', '\beta', '\gamma')

