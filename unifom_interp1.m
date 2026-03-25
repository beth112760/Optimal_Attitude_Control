format long
clear all; close all; clc;
%% Code to solve the attitude control problem
% Written by Bethany Hintz - 2 July 2025

% Initial Values
epsilon = 1;
N = 20;% Number of nodes
M = 3; % Number of elements
tf =  4.25; % Final time
ns = 1;
% Initial Parameteres
Neq = 14;
phi = pi;

load("CRBF_data.mat")

T1 = te(:);

T = te(:);

q_1 = Q1(:);
q_2 = Q2(:);
q_3 = Q3(:);
q_4 = Q4(:);
w_1 = W1(:);
w_2 = W2(:);
w_3 = W3(:);
u_1 = L5(:);
u_2 = L6(:);
u_3 = L7(:);

load("GPOPS.mat")
timevec = [solution.phase(1).time];


  
    q1vecc1 = [solution.phase(1).state(:,[1])];
 q2vecc1= [solution.phase(1).state(:,[2])];
q3vecc1 = [solution.phase(1).state(:,[3])];
 q4vecc1= [solution.phase(1).state(:,[4])];
 w1vecc1 = [solution.phase(1).state(:,[5])];
 w2vecc1= [solution.phase(1).state(:,[6])];
  w3vecc1 = [solution.phase(1).state(:,[7])];
  u1vecc1 = [solution.phase(1).control(:,[1])];
  u2vecc1 = [solution.phase(1).control(:,[2])];
  u3vecc1 = [solution.phase(1).control(:,[3])];


 q1vecc2 = [solution.phase(2).state(:,[1])];
 q2vecc2= [solution.phase(2).state(:,[2])];
 q3vecc2 = [solution.phase(2).state(:,[3])];
 q4vecc2= [solution.phase(2).state(:,[4])];
  w1vecc2 = [solution.phase(2).state(:,[5])];
  w2vecc2= [solution.phase(2).state(:,[6])];
  w3vecc2 = [solution.phase(2).state(:,[7])];
  u1vecc2 = [solution.phase(2).control(:,[1])];
  u2vecc2 = [solution.phase(2).control(:,[2])];
  u3vecc2 = [solution.phase(2).control(:,[3])];

      q1vecc3 = [solution.phase(3).state(:,[1])];
 q2vecc3= [solution.phase(3).state(:,[2])];
q3vecc3 = [solution.phase(3).state(:,[3])];
 q4vecc3= [solution.phase(3).state(:,[4])];
 w1vecc3 = [solution.phase(3).state(:,[5])];
 w2vecc3= [solution.phase(3).state(:,[6])];
  w3vecc3 = [solution.phase(3).state(:,[7])];
  u1vecc3 = [solution.phase(3).control(:,[1])];
  u2vecc3 = [solution.phase(3).control(:,[2])];
  u3vecc3 = [solution.phase(3).control(:,[3])];

  q1vecc = cat(1,q1vecc1,q1vecc2,q1vecc3);
 q2vecc= cat(1,q2vecc1,q2vecc2,q2vecc3);
q3vecc = cat(1,q3vecc1,q3vecc2,q3vecc3);
 q4vecc= cat(1,q4vecc1,q4vecc2,q4vecc3);
 w1vecc = cat(1,w1vecc1,w1vecc2,w1vecc3);
 w2vecc= cat(1,w2vecc1,w2vecc2,w2vecc3);
  w3vecc = cat(1,w3vecc1,w3vecc2,w3vecc3);
  u1vecc = cat(1,u1vecc1,u1vecc2,u1vecc3);
  u2vecc = cat(1,u2vecc1,u2vecc2,u2vecc3);
  u3vecc = cat(1,u3vecc1,u3vecc2,u3vecc3);
  timevec = [solution.phase(1).time;solution.phase(2).time;solution.phase(3).time];

  [~, ind] = unique(timevec);
duplicate_ind = flip(setdiff(1:size(timevec), ind));
duplicate_value = timevec(duplicate_ind);
for i = 1 : length(duplicate_ind)
    timevec(duplicate_ind(i))  = [];
    q1vecc(duplicate_ind(i)) = [];
    q2vecc(duplicate_ind(i)) = [];
    q3vecc(duplicate_ind(i)) = [];
    q4vecc(duplicate_ind(i)) = [];
    w1vecc(duplicate_ind(i)) = [];
    w2vecc(duplicate_ind(i)) = [];
    w3vecc(duplicate_ind(i)) = [];
    u1vecc(duplicate_ind(i)) = [];
    u2vecc(duplicate_ind(i)) = [];
    u3vecc(duplicate_ind(i)) = [];

end

q1vec = interp1(timevec,q1vecc,T,'linear');
q2vec = interp1(timevec,q2vecc,T,'linear');
q3vec = interp1(timevec,q3vecc,T,'linear');
q4vec =interp1(timevec,q4vecc,T,'linear');
w1vec = interp1(timevec,w1vecc,T,'linear');
w2vec = interp1(timevec,w2vecc,T,'linear');
w3vec = interp1(timevec,w3vecc,T,'linear');
u1vec = interp1(timevec,u1vecc,T,'linear');
u2vec = interp1(timevec,u2vecc,T,'linear');
u3vec = interp1(timevec,u3vecc,T,'linear');

abserror_q1 = (q_1 - q1vec);
abserror_q2 = (q_2 - q2vec);
abserror_q3 = (q_3 - q3vec);
abserror_q4 = (q_4 - q4vec);
abserror_w1 = (w_1 - w1vec);
abserror_w2 = (w_2 - w2vec);
abserror_w3 = (w_3 - w3vec);
abserror_u1 = (u_1 - u1vec);
abserror_u2 = (u_2 - u2vec);
abserror_u3 = (u_3 - u3vec);


%% Plot Results %%
figure(1)
plot(T,q_1,'*','LineWidth',1.1)
hold on
plot(T,q_2,'^','LineWidth',1.1)
plot(T,q_3,'*','LineWidth',1.5)
plot(T,q_4,'*','LineWidth',1.5)
plot(T,q1vec,'-','LineWidth',1.1)
plot(T,q2vec,'-','LineWidth',1.1)
plot(T,q3vec,'-','LineWidth',1.5)
plot(T,q4vec,'-','LineWidth',1.5)
xlabel('Time')
ylabel('Quarternion components')
title('Attitude Control - Quarternions')
legend('q1 (CRBF)','q2 (CRBF)','q3 (CRBF)','q4 (CRBF)','q1 (GPOPS-II)', 'q2 (GPOPS-II)', 'q3 (GPOPS-II)', 'q4 (GPOPS-II)')

figure(2)
plot(T,w_1,'*','LineWidth',1.1)
hold on
plot(T,w_2,'^','LineWidth',1.1)
plot(T,w_3,'*','LineWidth',1.5)
plot(T,w1vec,'-','LineWidth',1.1)
plot(T,w2vec,'-','LineWidth',1.1)
plot(T,w3vec,'-','LineWidth',1.5)
xlabel('Time')
ylabel('Angular Velocity')
title('Attitude Control - Angular Velocities')
legend('\omega_1 (CRBF)','\omega_2 (CRBF)','\omega_3 (CRBF)','\omega_1 (GPOPS-II)', '\omega_2 (GPOPS-II)', '\omega_3 (GPOPS-II)')

%Orientation angles
alpha = acos(1-2*q_2.^2-2*q_3.^2);
beta = acos(1-2*q_3.^2-2*q_1.^2);
gamma = acos(1-2*q_1.^2-2*q_2.^2);
alpha = rad2deg(alpha);
beta = rad2deg(beta);
gamma = rad2deg(gamma);

alphaa = acos(1-2*q2vec.^2-2*q3vec.^2);
betaa = acos(1-2*q3vec.^2-2*q1vec.^2);
gammaa = acos(1-2*q1vec.^2-2*q2vec.^2);
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
legend('\alpha (CRBF)','\beta (CRBF)','\gamma (CRBF)', '\alpha (GPOPS-II)','\beta (GPOPS-II)','\gamma (GPOPS-II)')

figure(4)
plot(T,u_1,'*','LineWidth',0.5)
hold on
plot(T,u_2,'^','LineWidth',1.5)
plot(T,u_3,'*','LineWidth',1.5)
plot(T,u1vec,'-','LineWidth',0.5)
plot(T,u2vec,'-','LineWidth',1.5)
plot(T,u3vec,'-','LineWidth',1.5)
xlabel('Time')
ylabel('Control Components along Body Axis')
title('Attitude Control - Control Components')
legend('u1 (CRBF)','u2 (CRBF)','u3 (CRBF)', 'u1 (GPOPS-II)','u2 (GPOPS-II)','u3 (GPOPS-II)')

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
plot(T,abs(abserror_u1),'-*','LineWidth',1.1)
hold on
plot(T,abs(abserror_u2),'-^','LineWidth',1.1)
plot(T,abs(abserror_u3),'-*','LineWidth',1.1)
xlabel('Time')
ylabel('Absolute Error')
title('Control Component Absolute Error')
legend('u_1', 'u_2', 'u_3')

figure;
plot(T,abs(abserror_alpha),'-*','LineWidth',1.1)
hold on
plot(T,abs(abserror_beta),'-^','LineWidth',1.1)
plot(T,abs(abserror_gamma),'-*','LineWidth',1.1)
xlabel('Time')
ylabel('Absolute Error')
title('Orientation Angles Absolute Error')
legend('\alpha', '\beta', '\gamma')

