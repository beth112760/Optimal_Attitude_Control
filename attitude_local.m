
function [R,J] = attitude_local(x,BC,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype)
global u1 u2 u3
%% State, Costates and BCs
x = reshape(x,[14*N,M]);
%ns = 1;
q1 = x(1:N,:);
q2 = x(N+1:2*N,:);
q3 = x(2*N+1:3*N,:);
q4 = x(3*N+1:4*N,:);
w1 = x(4*N+1:5*N,:);
w2 = x(5*N+1:6*N,:);
w3 = x(6*N+1:7*N,:);
% L1 = Lq1, L2 = Lq2, L3 = Lq3, L4 = Lq4, L5 = Lw1, L6 = Lw2, L7 = Lw3
L1 = x(7*N+1:8*N,:);
L2 = x(8*N+1:9*N,:);
L3 = x(9*N+1:10*N,:);
L4 = x(10*N+1:11*N,:);
L5 = x(11*N+1:12*N,:);
L6 = x(12*N+1:13*N,:);
L7 = x(13*N+1:14*N,:);

q10 = BC(1);
q20 = BC(2);
q30 = BC(3);
q40 = BC(4);
w10 = BC(5);
w20 = BC(6);
w30 = BC(7);
q1f = BC(8);
q2f = BC(9);
q3f = BC(10);
q4f = BC(11);
w1f = BC(12);
w2f = BC(13);
w3f = BC(14);



%% System of Equations
for k = 1 : M % for each cover

        %Residuals    
        R(1:N,k)       = D(:,:,k)*q1(:,k) - 0.5*(w1(:,k).*q4(:,k)-w2(:,k).*q3(:,k)+w3(:,k).*q2(:,k)); 
        R(N+1:2*N,k)   = D(:,:,k)*q2(:,k) - 0.5*(w1(:,k).*q3(:,k)+w2(:,k).*q4(:,k)-w3(:,k).*q1(:,k));
        R(2*N+1:3*N,k) = D(:,:,k)*q3(:,k) - 0.5*(-w1(:,k).*q2(:,k)+w2(:,k).*q1(:,k)+w3(:,k).*q4(:,k));
        R(3*N+1:4*N,k) = D(:,:,k)*q4(:,k) - 0.5*(-w1(:,k).*q1(:,k)-w2(:,k).*q2(:,k)-w3(:,k).*q3(:,k));    
        R(4*N+1:5*N,k) = D(:,:,k)*w1(:,k) - L5(:,k);
        R(5*N+1:6*N,k) = D(:,:,k)*w2(:,k) - L6(:,k);
        R(6*N+1:7*N,k) = D(:,:,k)*w3(:,k) - L7(:,k);
        % L1 = Lq1, L2 = Lq2, L3 = Lq3, L4 = Lq4, L5 = Lw1, L6 = Lw2, L7 = Lw3
        R(7*N+1:8*N,k) = D(:,:,k)*L1(:,k) + 0.5*(-L2(:,k).*w3(:,k) + L3(:,k).*w2(:,k) - L4(:,k).*w1(:,k));
        R(8*N+1:9*N,k) = D(:,:,k)*L2(:,k) + 0.5*(L1(:,k).*w3(:,k) - L4(:,k).*w2(:,k) - L3(:,k).*w1(:,k));
        R(9*N+1:10*N,k) = D(:,:,k)*L3(:,k) + 0.5*(-L4(:,k).*w3(:,k) - L1(:,k).*w2(:,k) + L2(:,k).*w1(:,k));
        R(10*N+1:11*N,k) = D(:,:,k)*L4(:,k) + 0.5*(L3(:,k).*w3(:,k) + L2(:,k).*w2(:,k) + L1(:,k).*w1(:,k));
        R(11*N+1:12*N,k) = D(:,:,k)*L5(:,k) + 0.5*(L1(:,k).*q4(:,k) + L2(:,k).*q3(:,k) - L3(:,k).*q2(:,k) - L4(:,k).*q1(:,k));
        R(12*N+1:13*N,k) = D(:,:,k)*L6(:,k) + 0.5*(L2(:,k).*q4(:,k) - L1(:,k).*q3(:,k) - L4(:,k).*q2(:,k) + L3(:,k).*q1(:,k));
        R(13*N+1:14*N,k) = D(:,:,k)*L7(:,k) + 0.5*(L3(:,k).*q4(:,k) - L4(:,k).*q3(:,k) + L1(:,k).*q2(:,k) - L2(:,k).*q1(:,k));

        % R(7*N+1:8*N,k) = D(:,:,k)*L1(:,k) + 0.5*(-w3(:,k) + w2(:,k) - w1(:,k));
        % R(8*N+1:9*N,k) = D(:,:,k)*L2(:,k) + 0.5*(w3(:,k) - w2(:,k) - w1(:,k));
        % R(9*N+1:10*N,k) = D(:,:,k)*L3(:,k) + 0.5*(-w3(:,k) - w2(:,k) + w1(:,k));
        % R(10*N+1:11*N,k) = D(:,:,k)*L4(:,k) + 0.5*(w3(:,k) + w2(:,k) + w1(:,k));
        % R(11*N+1:12*N,k) = D(:,:,k)*L5(:,k) + 0.5*(q4(:,k) + q3(:,k) - q2(:,k) - q1(:,k));
        % R(12*N+1:13*N,k) = D(:,:,k)*L6(:,k) + 0.5*(q4(:,k) - q3(:,k) - q2(:,k) + q1(:,k));
        % R(13*N+1:14*N,k) = D(:,:,k)*L7(:,k) + 0.5*(q4(:,k) - q3(:,k) + q2(:,k) - q1(:,k));

        %Jacobian
        dR1dq1 = D(:,:,k);      dR1dq2 = -0.5*w3(:,k).*eye(N);   dR1dq3 = 0.5.*w2(:,k).*eye(N);   dR1dq4 = -0.5.*w1(:,k).*eye(N);   dR1dw1 = -0.5.*q4(:,k).*eye(N);    dR1dw2 = 0.5.*q3(:,k).*eye(N);   dR1dw3 = -0.5.*q2(:,k).*eye(N);
            dR1dl1 = zeros(N);  dR1dl2 = zeros(N);       dR1dl3 = zeros(N);      dR1dl4 = zeros(N);       dR1dl5 = zeros(N);        dR1dl6 = zeros(N);      dR1dl7 = zeros(N);
        dR2dq1 = 0.5.*w3(:,k).*eye(N);   dR2dq2 = D(:,:,k);       dR2dq3 = -0.5.*w1(:,k).*eye(N);  dR2dq4 = -0.5.*w2(:,k).*eye(N);   dR2dw1 = -0.5.*q3(:,k).*eye(N);    dR2dw2 = -0.5.*q4(:,k).*eye(N);   dR2dw3 = 0.5.*q1(:,k).*eye(N);
            dR2dl1 = zeros(N);  dR2dl2 = zeros(N);       dR2dl3 = zeros(N);      dR2dl4 = zeros(N);       dR2dl5 = zeros(N);        dR2dl6 = zeros(N);      dR2dl7 = zeros(N);
        dR3dq1 = -0.5.*w2(:,k).*eye(N);  dR3dq2 = 0.5.*w1(:,k).*eye(N);    dR3dq3 = D(:,:,k);      dR3dq4 = -0.5.*w3(:,k).*eye(N);   dR3dw1 = 0.5.*q2(:,k).*eye(N);     dR3dw2 = -0.5.*q1(:,k).*eye(N);  dR3dw3 = -0.5.*q4(:,k).*eye(N);
            dR3dl1 = zeros(N);  dR3dl2 = zeros(N);       dR3dl3 = zeros(N);      dR3dl4 = zeros(N);       dR3dl5 = zeros(N);        dR3dl6 = zeros(N);      dR3dl7 = zeros(N);
        dR4dq1 = 0.5.*w1(:,k).*eye(N);   dR4dq2 = 0.5.*w2(:,k).*eye(N);    dR4dq3 = 0.5.*w3(:,k).*eye(N);   dR4dq4 = D(:,:,k);       dR4dw1 = 0.5.*q1(:,k).*eye(N);     dR4dw2 = 0.5.*q2(:,k).*eye(N);   dR4dw3 = 0.5.*q3(:,k).*eye(N);
            dR4dl1 = zeros(N);  dR4dl2 = zeros(N);       dR4dl3 = zeros(N);      dR4dl4 = zeros(N);       dR4dl5 = zeros(N);        dR4dl6 = zeros(N);      dR4dl7 = zeros(N);
        dR5dq1 = zeros(N);      dR5dq2 = zeros(N);       dR5dq3 = zeros(N);      dR5dq4 = zeros(N);       dR5dw1 = D(:,:,k);        dR5dw2 = zeros(N);      dR5dw3 = zeros(N);
            dR5dl1 = zeros(N);  dR5dl2 = zeros(N);       dR5dl3 = zeros(N);      dR5dl4 = zeros(N);       dR5dl5 = -eye(N);        dR5dl6 = zeros(N);      dR5dl7 = zeros(N);
        dR6dq1 = zeros(N);      dR6dq2 = zeros(N);       dR6dq3 = zeros(N);      dR6dq4 = zeros(N);       dR6dw1 = zeros(N);        dR6dw2 = D(:,:,k);      dR6dw3 = zeros(N);
            dR6dl1 = zeros(N);  dR6dl2 = zeros(N);       dR6dl3 = zeros(N);      dR6dl4 = zeros(N);       dR6dl5 = zeros(N);        dR6dl6 = -eye(N);      dR6dl7 = zeros(N);
        dR7dq1 = zeros(N);      dR7dq2 = zeros(N);       dR7dq3 = zeros(N);      dR7dq4 = zeros(N);       dR7dw1 = zeros(N);        dR7dw2 = zeros(N);      dR7dw3 = D(:,:,k);
            dR7dl1 = zeros(N);  dR7dl2 = zeros(N);       dR7dl3 = zeros(N);      dR7dl4 = zeros(N);       dR7dl5 = zeros(N);        dR7dl6 = zeros(N);      dR7dl7 = -eye(N);
        
        dR8dq1 = zeros(N);      dR8dq2 = zeros(N);       dR8dq3 = zeros(N);      dR8dq4 = zeros(N);       dR8dw1 = -0.5*L4(:,k).*eye(N);     dR8dw2 = 0.5*L3(:,k).*eye(N);    dR8dw3 = -0.5*L2(:,k).*eye(N);
            dR8dl1 = D(:,:,k);  dR8dl2 = -0.5.*w3(:,k).*eye(N);       dR8dl3 = 0.5.*w2(:,k).*eye(N);      dR8dl4 = -0.5.*w1(:,k).*eye(N);       dR8dl5 = zeros(N);        dR8dl6 = zeros(N);      dR8dl7 = zeros(N);
        %fixed R9 R11
       dR9dq1 = zeros(N);      dR9dq2 = zeros(N);       dR9dq3 = zeros(N);      dR9dq4 = zeros(N);       dR9dw1 = -0.5*L3(:,k).*eye(N);     dR9dw2 = -0.5*L4(:,k).*eye(N);   dR9dw3 = 0.5*L1(:,k).*eye(N);
            dR9dl1 = 0.5.*w3(:,k).*eye(N);  dR9dl2 = D(:,:,k);       dR9dl3 = -0.5.*w1(:,k).*eye(N);      dR9dl4 = -0.5.*w2(:,k).*eye(N);       dR9dl5 = zeros(N);        dR9dl6 = zeros(N);      dR9dl7 = zeros(N);
        dR10dq1 = zeros(N);     dR10dq2 = zeros(N);      dR10dq3 = zeros(N);     dR10dq4 = zeros(N);      dR10dw1 = 0.5*L2(:,k).*eye(N);     dR10dw2 = -0.5*L1(:,k).*eye(N);  dR10dw3 = -0.5*L4(:,k).*eye(N);
            dR10dl1 = -0.5.*w2(:,k).*eye(N); dR10dl2 = 0.5.*w1(:,k).*eye(N);      dR10dl3 = D(:,:,k);     dR10dl4 = -0.5.*w3(:,k).*eye(N);      dR10dl5 = zeros(N);       dR10dl6 = zeros(N);     dR10dl7 = zeros(N);
        dR11dq1 = zeros(N);     dR11dq2 = zeros(N);      dR11dq3 = zeros(N);     dR11dq4 = zeros(N);      dR11dw1 = 0.5*L1(:,k).*eye(N);     dR11dw2 = 0.5*L2(:,k).*eye(N);   dR11dw3 = 0.5*L3(:,k).*eye(N);
            dR11dl1 = 0.5.*w1(:,k).*eye(N); dR11dl2 = 0.5.*w2(:,k).*eye(N);      dR11dl3 = 0.5.*w3(:,k).*eye(N);     dR11dl4 = D(:,:,k);      dR11dl5 = zeros(N);       dR11dl6 = zeros(N);     dR11dl7 = zeros(N);
        dR12dq1 = -0.5*L4(:,k).*eye(N);  dR12dq2 = -0.5*L3(:,k).*eye(N);   dR12dq3 = 0.5*L2(:,k).*eye(N);   dR12dq4 = 0.5*L1(:,k).*eye(N);    dR12dw1 = zeros(N);       dR12dw2 = zeros(N);     dR12dw3 = zeros(N);
            dR12dl1 = 0.5.*q4(:,k).*eye(N); dR12dl2 = 0.5.*q3(:,k).*eye(N);      dR12dl3 = -0.5.*q2(:,k).*eye(N);     dR12dl4 = -0.5.*q1(:,k).*eye(N);      dR12dl5 = D(:,:,k);       dR12dl6 = zeros(N);     dR12dl7 = zeros(N);
        dR13dq1 = 0.5*L3(:,k).*eye(N);   dR13dq2 = -0.5*L4(:,k).*eye(N);   dR13dq3 = -0.5*L1(:,k).*eye(N);  dR13dq4 = 0.5*L2(:,k).*eye(N);    dR13dw1 = zeros(N);       dR13dw2 = zeros(N);     dR13dw3 = zeros(N);
            dR13dl1 = -0.5.*q3(:,k).*eye(N); dR13dl2 = 0.5.*q4(:,k).*eye(N);      dR13dl3 = 0.5.*q1(:,k).*eye(N);     dR13dl4 = -0.5.*q2(:,k).*eye(N);      dR13dl5 = zeros(N);       dR13dl6 = D(:,:,k);     dR13dl7 = zeros(N);
        dR14dq1 = -0.5*L2(:,k).*eye(N);  dR14dq2 = 0.5*L1(:,k).*eye(N);    dR14dq3 = -0.5*L4(:,k).*eye(N);  dR14dq4 = 0.5*L3(:,k).*eye(N);    dR14dw1 = zeros(N);       dR14dw2 = zeros(N);     dR14dw3 = zeros(N);
            dR14dl1 = 0.5.*q2(:,k).*eye(N); dR14dl2 = -0.5.*q1(:,k).*eye(N);      dR14dl3 = 0.5.*q4(:,k).*eye(N);     dR14dl4 = -0.5.*q3(:,k).*eye(N);      dR14dl5 = zeros(N);       dR14dl6 = zeros(N);     dR14dl7 = D(:,:,k);


  
    % Jacobian boundary conditions
    if k == 1
         dR1dq1(1,:) = [1,zeros(1,N-1)];  dR1dq2(1,:) = zeros(1,N);       dR1dq3(1,:) = zeros(1,N);       dR1dq4(1,:) = zeros(1,N);       dR1dw1(1,:) = zeros(1,N);       dR1dw2(1,:) = zeros(1,N);       dR1dw3(1,:) = zeros(1,N);... 
            dR1dl1(1,:) = zeros(1,N);     dR1dl2(1,:) = zeros(1,N);       dR1dl3(1,:) = zeros(1,N);       dR1dl4(1,:) = zeros(1,N);       dR1dl5(1,:) = zeros(1,N);       dR1dl6(1,:) = zeros(1,N);       dR1dl7(1,:) = zeros(1,N); 

        dR2dq1(1,:) = zeros(1,N);         dR2dq2(1,:) = [1,zeros(1,N-1)]; dR2dq3(1,:) = zeros(1,N);       dR2dq4(1,:) = zeros(1,N);       dR2dw1(1,:) = zeros(1,N);       dR2dw2(1,:) = zeros(1,N);       dR2dw3(1,:) = zeros(1,N);...
            dR2dl1(1,:) = zeros(1,N);     dR2dl2(1,:) = zeros(1,N);       dR2dl3(1,:) = zeros(1,N);       dR2dl4(1,:) = zeros(1,N);       dR2dl5(1,:) = zeros(1,N);       dR2dl6(1,:) = zeros(1,N);       dR2dl7(1,:) = zeros(1,N);

        dR3dq1(1,:) = zeros(1,N);         dR3dq2(1,:) = zeros(1,N);       dR3dq3(1,:) = [1,zeros(1,N-1)]; dR3dq4(1,:) = zeros(1,N);       dR3dw1(1,:) = zeros(1,N);       dR3dw2(1,:) = zeros(1,N);       dR3dw3(1,:) = zeros(1,N);...
            dR3dl1(1,:) = zeros(1,N);     dR3dl2(1,:) = zeros(1,N);       dR3dl3(1,:) = zeros(1,N);       dR3dl4(1,:) = zeros(1,N);       dR3dl5(1,:) = zeros(1,N);       dR3dl6(1,:) = zeros(1,N);       dR3dl7(1,:) = zeros(1,N);

        dR4dq1(1,:) = zeros(1,N);         dR4dq2(1,:) = zeros(1,N);       dR4dq3(1,:) = zeros(1,N);       dR4dq4(1,:) = [1,zeros(1,N-1)]; dR4dw1(1,:) = zeros(1,N);       dR4dw2(1,:) = zeros(1,N);       dR4dw3(1,:) = zeros(1,N);...
            dR4dl1(1,:) = zeros(1,N);     dR4dl2(1,:) = zeros(1,N);       dR4dl3(1,:) = zeros(1,N);       dR4dl4(1,:) = zeros(1,N);       dR4dl5(1,:) = zeros(1,N);       dR4dl6(1,:) = zeros(1,N);       dR4dl7(1,:) = zeros(1,N);

        dR5dq1(1,:) = zeros(1,N);         dR5dq2(1,:) = zeros(1,N);       dR5dq3(1,:) = zeros(1,N);       dR5dq4(1,:) = zeros(1,N);       dR5dw1(1,:) = [1,zeros(1,N-1)]; dR5dw2(1,:) = zeros(1,N);       dR5dw3(1,:) = zeros(1,N);...
            dR5dl1(1,:) = zeros(1,N);     dR5dl2(1,:) = zeros(1,N);       dR5dl3(1,:) = zeros(1,N);       dR5dl4(1,:) = zeros(1,N);       dR5dl5(1,:) = zeros(1,N);       dR5dl6(1,:) = zeros(1,N);       dR5dl7(1,:) = zeros(1,N);

        dR6dq1(1,:) = zeros(1,N);         dR6dq2(1,:) = zeros(1,N);       dR6dq3(1,:) = zeros(1,N);      dR6dq4(1,:) = zeros(1,N);        dR6dw1(1,:) = zeros(1,N);       dR6dw2(1,:) = [1,zeros(1,N-1)]; dR6dw3(1,:) = zeros(1,N); ...
            dR6dl1(1,:) = zeros(1,N);     dR6dl2(1,:) = zeros(1,N);       dR6dl3(1,:) = zeros(1,N);      dR6dl4(1,:) = zeros(1,N);        dR6dl5(1,:) = zeros(1,N);       dR6dl6(1,:) = zeros(1,N);       dR6dl7(1,:) = zeros(1,N);
    
        dR7dq1(1,:) = zeros(1,N);         dR7dq2(1,:) = zeros(1,N);       dR7dq3(1,:) = zeros(1,N);      dR7dq4(1,:) = zeros(1,N);        dR7dw1(1,:) = zeros(1,N);       dR7dw2(1,:) = zeros(1,N);       dR7dw3(1,:) = [1,zeros(1,N-1)]; ...
            dR7dl1(1,:) = zeros(1,N);     dR7dl2(1,:) = zeros(1,N);       dR7dl3(1,:) = zeros(1,N);      dR7dl4(1,:) = zeros(1,N);        dR7dl5(1,:) = zeros(1,N);       dR7dl6(1,:) = zeros(1,N);       dR7dl7(1,:) = zeros(1,N);

    end
    if k == M
        %         dR7dx1(N,:) = [zeros(1,N-1),1]; dR7dx2(N,:) = zeros(1,N);       dR7dy1(N,:) = zeros(1,N); dR7dy2(N,:) = zeros(1,N); dR7dz1(N,:) = zeros(1,N); dR7dz2(N,:) = zeros(1,N); ... 
        %     dR7dl1(N,:) = zeros(1,N); dR7dl2(N,:) = zeros(1,N);dR7dl3(N,:) = zeros(1,N);dR7dl4(N,:) = zeros(1,N);dR7dl5(N,:) = zeros(1,N);dR7dl6(N,:) = zeros(1,N); 
        % 
        % dR8dx1(N,:) = zeros(1,N);       dR8dx2(N,:) = [zeros(1,N-1),1]; dR8dy1(N,:) = zeros(1,N); dR8dy2(N,:) = zeros(1,N); dR8dz1(N,:) = zeros(1,N); dR8dz2(N,:) = zeros(1,N);...
        %     dR8dl1(N,:) = zeros(1,N); dR8dl2(N,:) = zeros(1,N);dR8dl3(N,:) = zeros(1,N);dR8dl4(N,:) = zeros(1,N);dR8dl5(N,:) = zeros(1,N);dR8dl6(N,:) = zeros(1,N);
        % 
        % dR9dx1(N,:) = zeros(1,N); dR9dx2(N,:) = zeros(1,N);       dR9dy1(N,:) = [zeros(1,N-1),1]; dR9dy2(N,:) = zeros(1,N); dR9dz1(N,:) = zeros(1,N); dR9dz2(N,:) = zeros(1,N);...
        %     dR9dl1(N,:) = zeros(1,N); dR9dl2(N,:) = zeros(1,N);dR9dl3(N,:) = zeros(1,N);dR9dl4(N,:) = zeros(1,N);dR9dl5(N,:) = zeros(1,N);dR9dl6(N,:) = zeros(1,N);
        % 
        % dR10dx1(N,:) = zeros(1,N);       dR10dx2(N,:) = zeros(1,N); dR10dy1(N,:) = zeros(1,N); dR10dy2(N,:) = [zeros(1,N-1),1]; dR10dz1(N,:) = zeros(1,N); dR10dz2(N,:) = zeros(1,N);...
        %     dR10dl1(N,:) = zeros(1,N); dR10dl2(N,:) = zeros(1,N);dR10dl3(N,:) = zeros(1,N);dR10dl4(N,:) = zeros(1,N);dR10dl5(N,:) = zeros(1,N);dR10dl6(N,:) = zeros(1,N);
        % 
        % dR11dx1(N,:) = zeros(1,N); dR11dx2(N,:) = zeros(1,N);       dR11dy1(N,:) = zeros(1,N); dR11dy2(N,:) = zeros(1,N); dR11dz1(N,:) = [zeros(1,N-1),1]; dR11dz2(N,:) = zeros(1,N);...
        %     dR11dl1(N,:) = zeros(1,N); dR11dl2(N,:) = zeros(1,N);dR11dl3(N,:) = zeros(1,N);dR11dl4(N,:) = zeros(1,N);dR11dl5(N,:) = zeros(1,N);dR11dl6(N,:) = zeros(1,N);
        % 
        % dR12dx1(N,:) = zeros(1,N);       dR12dx2(N,:) = zeros(1,N); dR12dy1(N,:) = zeros(1,N); dR12dy2(N,:) = zeros(1,N); dR12dz1(N,:) = zeros(1,N); dR12dz2(N,:) = [zeros(1,N-1),1]; ...
        %     dR12dl1(N,:) = zeros(1,N); dR12dl2(N,:) = zeros(1,N);dR12dl3(N,:) = zeros(1,N);dR12dl4(N,:) = zeros(1,N);dR12dl5(N,:) = zeros(1,N);dR12dl6(N,:) = zeros(1,N);
        % 
        % dR8dq1(1,:) = [1,zeros(1,N-1)]; dR8dq2(1,:) = zeros(1,N);       dR8dq3(1,:) = zeros(1,N); dR8dq4(1,:) = zeros(1,N); dR8dw1(1,:) = zeros(1,N); dR8dw2(1,:) = zeros(1,N); dR8dw3(1,:) = zeros(1,N);... 
        %     dR8dl1(1,:) = zeros(1,N); dR8dl2(1,:) = zeros(1,N);dR8dl3(1,:) = zeros(1,N);dR8dl4(1,:) = zeros(1,N);dR8dl5(1,:) = zeros(1,N);dR8dl6(1,:) = zeros(1,N); dR8dl7(1,:) = zeros(1,N);
        % 
        % dR9dq1(1,:) = zeros(1,N);       dR9dq2(1,:) = [1,zeros(1,N-1)]; dR9dq3(1,:) = zeros(1,N); dR9dq4(1,:) = zeros(1,N); dR9dw1(1,:) = zeros(1,N); dR9dw2(1,:) = zeros(1,N); dR9dw3(1,:) = zeros(1,N);...
        %     dR9dl1(1,:) = zeros(1,N); dR9dl2(1,:) = zeros(1,N);dR9dl3(1,:) = zeros(1,N);dR9dl4(1,:) = zeros(1,N);dR9dl5(1,:) = zeros(1,N);dR9dl6(1,:) = zeros(1,N); dR9dl7(1,:) = zeros(1,N);
        % 
        % dR10dq1(1,:) = zeros(1,N); dR10dq2(1,:) = zeros(1,N);       dR10dq3(1,:) = [1,zeros(1,N-1)]; dR10dq4(1,:) = zeros(1,N); dR10dw1(1,:) = zeros(1,N); dR10dw2(1,:) = zeros(1,N); dR10dw3(1,:) = zeros(1,N);...
        %     dR10dl1(1,:) = zeros(1,N); dR10dl2(1,:) = zeros(1,N);dR10dl3(1,:) = zeros(1,N);dR10dl4(1,:) = zeros(1,N);dR10dl5(1,:) = zeros(1,N);dR10dl6(1,:) = zeros(1,N); dR10dl7(1,:) = zeros(1,N);
        % 
        % dR11dq1(1,:) = zeros(1,N);       dR11dq2(1,:) = zeros(1,N); dR11dq3(1,:) = zeros(1,N); dR11dq4(1,:) = [1,zeros(1,N-1)]; dR11dw1(1,:) = zeros(1,N); dR11dw2(1,:) = zeros(1,N); dR11dw3(1,:) = zeros(1,N);...
        %     dR11dl1(1,:) = zeros(1,N); dR11dl2(1,:) = zeros(1,N);dR11dl3(1,:) = zeros(1,N);dR11dl4(1,:) = zeros(1,N);dR11dl5(1,:) = zeros(1,N);dR11dl6(1,:) = zeros(1,N); dR11dl7(1,:) = zeros(1,N);
        % 
        % dR12dq1(1,:) = zeros(1,N); dR12dq2(1,:) = zeros(1,N);       dR12dq3(1,:) = zeros(1,N); dR12dq4(1,:) = zeros(1,N); dR12dw1(1,:) = [1,zeros(1,N-1)]; dR12dw2(1,:) = zeros(1,N); dR12dw3(1,:) = zeros(1,N);...
        %     dR12dl1(1,:) = zeros(1,N); dR12dl2(1,:) = zeros(1,N);dR12dl3(1,:) = zeros(1,N);dR12dl4(1,:) = zeros(1,N);dR12dl5(1,:) = zeros(1,N);dR12dl6(1,:) = zeros(1,N); dR12dl7(1,:) = zeros(1,N);
        % 
        % dR13dq1(1,:) = zeros(1,N);       dR13dq2(1,:) = zeros(1,N); dR13dq3(1,:) = zeros(1,N); dR13dq4(1,:) = zeros(1,N); dR13dw1(1,:) = zeros(1,N); dR13dw2(1,:) = [1,zeros(1,N-1)]; dR13dw3(1,:) = zeros(1,N); ...
        %     dR13dl1(1,:) = zeros(1,N); dR13dl2(1,:) = zeros(1,N);dR13dl3(1,:) = zeros(1,N);dR13dl4(1,:) = zeros(1,N);dR13dl5(1,:) = zeros(1,N);dR13dl6(1,:) = zeros(1,N); dR13dl7(1,:) = zeros(1,N);
        % dR14dq1(1,:) = zeros(1,N);       dR14dq2(1,:) = zeros(1,N); dR14dq3(1,:) = zeros(1,N); dR14dq4(1,:) = zeros(1,N); dR14dw1(1,:) = zeros(1,N); dR14dw2(1,:) = zeros(1,N); dR14dw3(1,:) = [1,zeros(1,N-1)]; ...
        %     dR14dl1(1,:) = zeros(1,N); dR14dl2(1,:) = zeros(1,N);dR14dl3(1,:) = zeros(1,N);dR14dl4(1,:) = zeros(1,N);dR14dl5(1,:) = zeros(1,N);dR14dl6(1,:) = zeros(1,N); dR14dl7(1,:) = zeros(1,N);
        % 



        dR7dq1(1,:) = [1,zeros(1,N-1)]; dR7dq2(1,:) = zeros(1,N);       dR7dq3(1,:) = zeros(1,N); dR7dq4(1,:) = zeros(1,N); dR7dw1(1,:) = zeros(1,N); dR7dw2(1,:) = zeros(1,N); dR7dw3(1,:) = zeros(1,N);... 
            dR7dl1(1,:) = zeros(1,N); dR7dl2(1,:) = zeros(1,N);dR7dl3(1,:) = zeros(1,N);dR7dl4(1,:) = zeros(1,N);dR7dl5(1,:) = zeros(1,N);dR7dl6(1,:) = zeros(1,N); dR7dl7(1,:) = zeros(1,N);

        dR8dq1(1,:) = zeros(1,N);       dR8dq2(1,:) = [1,zeros(1,N-1)]; dR8dq3(1,:) = zeros(1,N); dR8dq4(1,:) = zeros(1,N); dR8dw1(1,:) = zeros(1,N); dR8dw2(1,:) = zeros(1,N); dR8dw3(1,:) = zeros(1,N);...
            dR8dl1(1,:) = zeros(1,N); dR8dl2(1,:) = zeros(1,N);dR8dl3(1,:) = zeros(1,N);dR8dl4(1,:) = zeros(1,N);dR8dl5(1,:) = zeros(1,N);dR8dl6(1,:) = zeros(1,N); dR8dl7(1,:) = zeros(1,N);

        dR9dq1(1,:) = zeros(1,N); dR9dq2(1,:) = zeros(1,N);       dR9dq3(1,:) = [1,zeros(1,N-1)]; dR9dq4(1,:) = zeros(1,N); dR9dw1(1,:) = zeros(1,N); dR9dw2(1,:) = zeros(1,N); dR9dw3(1,:) = zeros(1,N);...
            dR9dl1(1,:) = zeros(1,N); dR9dl2(1,:) = zeros(1,N);dR9dl3(1,:) = zeros(1,N);dR9dl4(1,:) = zeros(1,N);dR9dl5(1,:) = zeros(1,N);dR9dl6(1,:) = zeros(1,N); dR9dl7(1,:) = zeros(1,N);

        dR10dq1(1,:) = zeros(1,N);       dR10dq2(1,:) = zeros(1,N); dR10dq3(1,:) = zeros(1,N); dR10dq4(1,:) = [1,zeros(1,N-1)]; dR10dw1(1,:) = zeros(1,N); dR10dw2(1,:) = zeros(1,N); dR10dw3(1,:) = zeros(1,N);...
            dR10dl1(1,:) = zeros(1,N); dR10dl2(1,:) = zeros(1,N);dR10dl3(1,:) = zeros(1,N);dR10dl4(1,:) = zeros(1,N);dR10dl5(1,:) = zeros(1,N);dR10dl6(1,:) = zeros(1,N); dR10dl7(1,:) = zeros(1,N);

        dR11dq1(1,:) = zeros(1,N); dR11dq2(1,:) = zeros(1,N);       dR11dq3(1,:) = zeros(1,N); dR11dq4(1,:) = zeros(1,N); dR11dw1(1,:) = [1,zeros(1,N-1)]; dR11dw2(1,:) = zeros(1,N); dR11dw3(1,:) = zeros(1,N);...
            dR11dl1(1,:) = zeros(1,N); dR11dl2(1,:) = zeros(1,N);dR11dl3(1,:) = zeros(1,N);dR11dl4(1,:) = zeros(1,N);dR11dl5(1,:) = zeros(1,N);dR11dl6(1,:) = zeros(1,N); dR11dl7(1,:) = zeros(1,N);

        dR12dq1(1,:) = zeros(1,N);       dR12dq2(1,:) = zeros(1,N); dR12dq3(1,:) = zeros(1,N); dR12dq4(1,:) = zeros(1,N); dR12dw1(1,:) = zeros(1,N); dR12dw2(1,:) = [1,zeros(1,N-1)]; dR12dw3(1,:) = zeros(1,N); ...
            dR12dl1(1,:) = zeros(1,N); dR12dl2(1,:) = zeros(1,N);dR12dl3(1,:) = zeros(1,N);dR12dl4(1,:) = zeros(1,N);dR12dl5(1,:) = zeros(1,N);dR12dl6(1,:) = zeros(1,N); dR12dl7(1,:) = zeros(1,N);
        dR13dq1(1,:) = zeros(1,N);       dR13dq2(1,:) = zeros(1,N); dR13dq3(1,:) = zeros(1,N); dR13dq4(1,:) = zeros(1,N); dR13dw1(1,:) = zeros(1,N); dR13dw2(1,:) = zeros(1,N); dR13dw3(1,:) = [1,zeros(1,N-1)]; ...
            dR13dl1(1,:) = zeros(1,N); dR13dl2(1,:) = zeros(1,N);dR13dl3(1,:) = zeros(1,N);dR13dl4(1,:) = zeros(1,N);dR13dl5(1,:) = zeros(1,N);dR13dl6(1,:) = zeros(1,N); dR13dl7(1,:) = zeros(1,N);

     end
    J_local(:,:,k) = [dR1dq1, dR1dq2, dR1dq3, dR1dq4 dR1dw1, dR1dw2, dR1dw3, dR1dl1, dR1dl2, dR1dl3, dR1dl4, dR1dl5, dR1dl6, dR1dl7;...
                      dR2dq1, dR2dq2, dR2dq3, dR2dq4 dR2dw1, dR2dw2, dR2dw3, dR2dl1, dR2dl2, dR2dl3, dR2dl4, dR2dl5, dR2dl6, dR2dl7;...
                      dR3dq1, dR3dq2, dR3dq3, dR3dq4 dR3dw1, dR3dw2, dR3dw3, dR3dl1, dR3dl2, dR3dl3, dR3dl4, dR3dl5, dR3dl6, dR3dl7;...
                      dR4dq1, dR4dq2, dR4dq3, dR4dq4 dR4dw1, dR4dw2, dR4dw3, dR4dl1, dR4dl2, dR4dl3, dR4dl4, dR4dl5, dR4dl6, dR4dl7;...
                      dR5dq1, dR5dq2, dR5dq3, dR5dq4 dR5dw1, dR5dw2, dR5dw3, dR5dl1, dR5dl2, dR5dl3, dR5dl4, dR5dl5, dR5dl6, dR5dl7;...
                      dR6dq1, dR6dq2, dR6dq3, dR6dq4 dR6dw1, dR6dw2, dR6dw3, dR6dl1, dR6dl2, dR6dl3, dR6dl4, dR6dl5, dR6dl6, dR6dl7;...
                      dR7dq1, dR7dq2, dR7dq3, dR7dq4 dR7dw1, dR7dw2, dR7dw3, dR7dl1, dR7dl2, dR7dl3, dR7dl4, dR7dl5, dR7dl6, dR7dl7;...
                      dR8dq1, dR8dq2, dR8dq3, dR8dq4 dR8dw1, dR8dw2, dR8dw3, dR8dl1, dR8dl2, dR8dl3, dR8dl4, dR8dl5, dR8dl6, dR8dl7;...
                      dR9dq1, dR9dq2, dR9dq3, dR9dq4 dR9dw1, dR9dw2, dR9dw3, dR9dl1, dR9dl2, dR9dl3, dR9dl4, dR9dl5, dR9dl6, dR9dl7;...
                      dR10dq1, dR10dq2, dR10dq3, dR10dq4 dR10dw1, dR10dw2, dR10dw3, dR10dl1, dR10dl2, dR10dl3, dR10dl4, dR10dl5, dR10dl6, dR10dl7;...
                      dR11dq1, dR11dq2, dR11dq3, dR11dq4 dR11dw1, dR11dw2, dR11dw3, dR11dl1, dR11dl2, dR11dl3, dR11dl4, dR11dl5, dR11dl6, dR11dl7;...
                      dR12dq1, dR12dq2, dR12dq3, dR12dq4 dR12dw1, dR12dw2, dR12dw3, dR12dl1, dR12dl2, dR12dl3, dR12dl4, dR12dl5, dR12dl6, dR12dl7; ...
                      dR13dq1, dR13dq2, dR13dq3, dR13dq4 dR13dw1, dR13dw2, dR13dw3, dR13dl1, dR13dl2, dR13dl3, dR13dl4, dR13dl5, dR13dl6, dR13dl7; ...
                      dR14dq1, dR14dq2, dR14dq3, dR14dq4 dR14dw1, dR14dw2, dR14dw3, dR14dl1, dR14dl2, dR14dl3, dR14dl4, dR14dl5, dR14dl6, dR14dl7];


    % Indices of states and segments
    Elb(k) = (k-1) * LJ_local+1;        %left boundary of each segment (Element left boundary)
    Erb(k) = Elb(k) + LJ_local- 1;      %right boundary of each segment
    for ii = 1 : Neq
        SLb(ii,:,k) = Elb(k)+ (ii-1)*N; %left boundary of each state ii for each element k (State left boundary)
        SRb(ii,:,k) = Elb(k) + ii*N-1;
        slb(ii) = (ii-1)*N+1;           %Left boundary of each state ii for the "local matrix"
        srb(ii) = ii*N;                 %Right boundary of each state ii for the "local matrix"
    end
    
    %Global Jacobian matrix
    J(Elb(k):Erb(k),Elb(k):Erb(k)) = J_local(:,:,k); %Global matrix

    % Continuity Conditions of the Jacobian Matrix
    if k > 1
        for ii = 1 : Neq
            J(SLb(ii,:,k):SLb(ii,:,k)+ns,:)                           = zeros(1+ns,length(J));
            J(SLb(ii,:,k):SLb(ii,:,k)+ns,SLb(ii,:,k-1):SRb(ii,:,k-1)) = [zeros(1,N-1),-1];
            J(SLb(ii,:,k):SLb(ii,:,k)+ns,SLb(ii,:,k):SRb(ii,:,k)) = [1,zeros(1,N-1)];
        end
    end

        if k > 1 
        R(1,k)     = q1(1,k) - q1(end,k-1);
        R(N+1,k)   = q2(1,k) - q2(end,k-1);
        R(2*N+1,k) = q3(1,k) - q3(end,k-1);
        R(3*N+1,k) = q4(1,k) - q4(end,k-1);
        R(4*N+1,k) = w1(1,k) - w1(end,k-1);
        R(5*N+1,k) = w2(1,k) - w2(end,k-1);
        R(6*N+1,k) = w3(1,k) - w3(end,k-1);
        R(7*N+1,k) = L1(1,k) - L1(end,k-1);
        R(8*N+1,k) = L2(1,k) - L2(end,k-1);
        R(9*N+1,k) = L3(1,k) - L3(end,k-1);
        R(10*N+1,k) = L4(1,k) - L4(end,k-1);
        R(11*N+1,k) = L5(1,k) - L5(end,k-1);
        R(12*N+1,k) = L6(1,k) - L6(end,k-1);
        R(13*N+1,k) = L7(1,k) - L7(end,k-1);
        end
end
%% Boundary conditions
switch BCtype
    case "fixed" %fixed final state
 
        R(1,1)   = q1(1,1) - q10;
        R(N+1,1) = q2(1,1) - q20;
        R(7*N+1,1)   = q1(end,end) - q1f;
        R(8*N+1,1)   = q2(end,end) - q2f;
        
        R(2*N+1,1)   = q3(1,1) - q30;
        R(3*N+1,1) = q4(1,1) - q40;
        R(9*N+1,1)   = q3(end,end) - q3f;
        R(10*N+1,1)   = q4(end,end) - q4f;

        R(4*N+1,1)   = w1(1,1) - w10;
        R(5*N+1,1) = w2(1,1) - w20;
        R(6*N+1,1) = w3(1,1) - w30;
        R(11*N+1,1)   = w1(end,end) - w1f;
        R(12*N+1,1)   = w2(end,end) - w2f;
        R(13*N+1,1)   = w3(end,end) - w3f;

    case "free" %free final state
end

%% Output Residual Vector
% Function output
 R = R(:);  %[all state segment1; all states segment2; ... ; all states last segment];

if ip == 1
    R = 1/2*(R'*R);
end


