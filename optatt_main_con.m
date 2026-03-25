format long
clear all; close all; clc;
%% Code to solve the attitude control problem
% Written by Bethany Hintz - 2 July 2025

% Initial Values
epsilon = 1;
N = 54;% Number of nodes
M = 3; % Number of elements
tf =  4.25; % Final time
ns = 1;
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
t1 = tf/2-0.209;
t2 = tf/2+0.209;
%% Approximation matrix
D = zeros(N,N,M);
phi = zeros(N,N,M);
phid = zeros(N,N,M);

if mod(ns, 2) == 0 %if ns is even

%Divide time domain to accomodate shared nodes
    for k = 1:M
        if k == 1
            t0e            = 0; % Initial time of each segment (element)
            tfe            = t1;     % Final time of each segment (element)
            temp(:,:,k)      = linspace(t0e,tfe,N-floor(ns/2));  
         elseif k == M
             t0e            = t2; % Initial time of each segment (element)
             tfe            = tf;     % Final time of each segment (element)
             temp2(:,:,k)      = linspace(t0e,tfe,N-(floor(ns/2))+1); 
        elseif k == 2 
            t0e            = t1; % Initial time of each segment (element)
            tfe            = (k-1)*(t2-t1)/(M-2)+t1; (k+1)*t2/M;     % Final time of each segment (element)
            temp1(:,:,k)      = linspace(t0e,tfe,N-floor(ns/2)-(floor(ns/2))+1);
             elseif k == M-1
            t0e            = (k-2)*(t2-t1)/(M-2)+t1; % Initial time of each segment (element)
            tfe            = t2;     % Final time of each segment (element)
            temp1(:,:,k)      = linspace(t0e,tfe,N-floor(ns/2)-(floor(ns/2))+1); 
          else
            t0e            = (k-2)*(t2-t1)/(M-2)+t1; % Initial time of each segment (element)
            tfe            = (k-1)*(t2-t1)/(M-2)+t1;     % Final time of each segment (element)
            temp1(:,:,k)      = linspace(t0e,tfe,N-floor(ns/2)-(floor(ns/2))+1);
       end
    end

    %Extract shared nodes
    for k = 1
        cat_temp(:,:,k) = cat(2,temp(:,:,k),temp1(:,2:ns-floor(ns/2)+1,k+1));
    end
    if M > 3
        for k = 2 % First "middle" element
            cat_temp1(:,:,k) = cat(2,temp(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp1(:,2:ns-floor(ns/2)+1,k+1));
        end
        if M == 4
            for k = 3
                cat_temp3(:,:,k) = cat(2,temp1(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp2(:,2:ns-floor(ns/2)+1,k+1));
            end
        else
            for k = 3:M-2 
                cat_temp3(:,:,k) = cat(2,temp1(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp1(:,2:ns-floor(ns/2)+1,k+1));
            end
        end
        for k = M-1 % Last "middle" element
                cat_temp3(:,:,k) = cat(2,temp1(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp2(:,2:ns-floor(ns/2)+1,k+1));        
        end
    else
        for k = 2
            cat_temp1(:,:,k) = cat(2,temp(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp2(:,2:ns-floor(ns/2)+1,k+1));
        end
    end
    for k = M
        cat_temp2(:,:,k) = cat(2,temp1(:,end-ns+floor(ns/2)+1:end-1,k-1),temp2(:,:,k));
    end

    %Permute and Reshape 3d arrays
    series1 = reshape(cat_temp,[length(cat_temp),1]);
    C = permute(cat_temp1,[1 3 2]);
    C = reshape(C,[],size(cat_temp1,2),1);
    series2 = nonzeros(C);
    C = permute(cat_temp2,[1 3 2]);
    C = reshape(C,[],size(cat_temp2,2),1);
    series3 = nonzeros(C);
    if M >3
    C = permute(cat_temp3,[1 3 2]);
    C = reshape(C,[],size(cat_temp3,2),1)';
    series4 = nonzeros(C);
    end

    %Build time vector with shared nodes
    if M > 3
        te = cat(1,series1,series2,series4,series3);
    else
        te=cat(1,series1,series2,series3);
    end
    te = reshape(te,1,N,M);
else %if ns is odd

    % Middle cover region
    tempy = linspace(t1, t2, (M-2)*(N-2*floor(ns/2)));

%Divide time domain to accomodate shared nodes
    for k = 1:M
        if k == 1 
            t0e            = 0; % Initial time of each segment (element)
            tfe            = t1;     % Final time of each segment (element)
            temp(:,:,k)      = linspace(t0e,tfe,N-floor(ns/2)); 
        elseif k == M
            t0e            = t2; % Initial time of each segment (element)
            tfe            = tf;     % Final time of each segment (element)
            temp(:,:,k)      = linspace(t0e,tfe,N-floor(ns/2));             
        
        elseif k == 2 
            t0e            = t1; % Initial time of each segment (element)
            tfe            = (k-1)*(t2-t1)/(M-2)+t1;     % Final time of each segment (element)
            temp1(:,:,k)      = linspace(t0e,tfe,N-2*floor(ns/2));       
        elseif k == M-1
            t0e            = (k-2)*(t2-t1)/(M-2)+t1; % Initial time of each segment (element)
            tfe            = t2;     % Final time of each segment (element)
            temp1(:,:,k)      = linspace(t0e,tfe,N-2*floor(ns/2)); 
             
        else
            t0e            = (k-2)*(t2-t1)/(M-2)+t1; % Initial time of each segment (element)
            tfe            = (k-1)*(t2-t1)/(M-2)+t1;     % Final time of each segment (element)
            temp1(:,:,k)      = linspace(t0e,tfe,N-2*floor(ns/2)); 
        end
    end

    %Extract shared nodes
    for k = 1 % Append Element #1 with the shared nodes existing in Element #2
        cat_temp(:,:,k) = cat(2,temp(:,:,k),temp1(:,2:ns-floor(ns/2),k+1));
    end
    % To account for the cases where there are multiple "middle" elements
    if M > 3
        for k = 2 % First "middle" element
            cat_temp1(:,:,k) = cat(2,temp(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp1(:,2:ns-floor(ns/2),k+1));
        end
        for k = 3:M-2 
            cat_temp1(:,:,k) = cat(2,temp1(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp1(:,2:ns-floor(ns/2),k+1));
        end
        for k = M-1 % Last "middle" element
                cat_temp1(:,:,k) = cat(2,temp1(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp(:,2:ns-floor(ns/2),k+1));        
        end
    % To account for the case where there is a single "middle" element (i.e. M = 3)    
    else
        for k = 2:M-1
            cat_temp1(:,:,k) = cat(2,temp(:,end-ns+floor(ns/2)+1:end-1,k-1),temp1(:,:,k),temp(:,2:ns-floor(ns/2),k+1));
        end
    end
    % Append last Element with shared nodes from preceding Element
    for k = M
        cat_temp2(:,:,k) = cat(2,temp1(:,end-ns+floor(ns/2)+1:end-1,k-1),temp(:,:,k));
    end
    
    %Permute and Reshape 3d arrays
    series1 = reshape(cat_temp,[length(cat_temp),1]);
    C = permute(cat_temp1,[1 3 2]);
    C = reshape(C,[],size(cat_temp1,2),1)';
    series2 = nonzeros(C);
    C = permute(cat_temp2,[1 3 2]);
    C = reshape(C,[],size(cat_temp2,2),1)';
    series3 = nonzeros(C);
    
    %Build time vector with shared nodes
    te = cat(1,series1,series2,series3);
    te = reshape(te,1,N,M);
end

for k = 1 : M
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end

    D(:,:,k) = phid(:,:,k)/(phi(:,:,k));

%% CRBF Solution using fsolve
ip=0;
Jlocal_pattern = ones(Neq*N);
J_pattern = kron(eye(M),Jlocal_pattern);

J_local = zeros(Neq*N,Neq*N,M);
J    = zeros(Neq*N*M, Neq*N*M); % Local matrices are to be concatenated sequentially

Elb   = zeros(1,M);
Erb   = Elb;
SLb   = zeros(Neq,1,M);

SRb   = SLb;
R = zeros(N,M);
sizeJ_local = size(J_local);
LJ_local    = sizeJ_local(1);

f = @(X) attitude_local3(X,BC,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
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

% %% Clean Data %%
T = te(:);

q1 = Q1(:);
q2 = Q2(:);
q3 = Q3(:);
q4 = Q4(:);
w1 = W1(:);
w2 = W2(:);
w3 = W3(:);
u1 = L5(:);
u2 = L6(:);
u3 = L7(:);

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
uu1 = LL5(:);
uu2 = LL6(:);
uu3 = LL7(:);

abserror_q1 = (q1 - qq1);
abserror_q2 = (q2 - qq2);
abserror_q3 = (q3 - qq3);
abserror_q4 = (q4 - qq4);
abserror_w1 = (w1 - ww1);
abserror_w2 = (w2 - ww2);
abserror_w3 = (w3 - ww3);
abserror_u1 = (u1 - uu1);
abserror_u2 = (u2 - uu2);
abserror_u3 = (u3 - uu3);


%% Plot Results %%
figure(1)
plot(T,q1,'*','LineWidth',1.1)
hold on
plot(T,q2,'^','LineWidth',1.1)
plot(T,q3,'*','LineWidth',1.5)
plot(T,q4,'*','LineWidth',1.5)
xlabel('Time')
ylabel('Quarternion components')
title('Attitude Control - Quarternions')
legend('q1 (CRBF)','q2 (CRBF)','q3 (CRBF)','q4 (CRBF)','q1 (Exact)', 'q2 (Exact)', 'q3 (Exact)', 'q4 (Exact)')

figure
plot(T,w1,'*','LineWidth',1.1)
hold on
plot(T,w2,'^','LineWidth',1.1)
plot(T,w3,'*','LineWidth',1.5)
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
xlabel('Time')
ylabel('Orientation angle (deg)')
title('Attitude Control - Orientation Angles')
legend('\alpha (CRBF)','\beta (CRBF)','\gamma (CRBF)', '\alpha (Exact)','\beta (Exact)','\gamma (Exact)')

figure(4)
plot(T,u1,'*','LineWidth',0.5)
hold on
plot(T,u2,'^','LineWidth',1.5)
plot(T,u3,'*','LineWidth',1.5)
xlabel('Time')
ylabel('Control Components along Body Axis')
title('Attitude Control - Control Components')
legend('u1 (CRBF)','u2 (CRBF)','u3 (CRBF)', 'u1 (Exact)','u2 (Exact)','u3 (Exact)')


