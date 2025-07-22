function states = attdynamics(x)

q1 = x(1);
q2 = x(2);
q3 = x(3);
q4 = x(4);
w1 = x(5);
w2 = x(6);
w3 = x(7);

lambda1 = x(8);
lambda2 = x(9);
lambda3 = x(10);
lambda4 = x(11);
lambda5 = x(12);
lambda6 = x(13);
lambda7 = x(14);

q1dot = 0.5*(w1*q4 - w2*q3 + w3*q2);
q2dot = 0.5*(w1*q3 + w2*q4 - w3*q1);
q3dot= 0.5*(-w1*q2 + w2*q1 + w3*q4);
q4dot = 0.5*(-w1*q1 - w2*q2 - w3*q3);
w1dot= lambda5;
w2dot = lambda6;
w3dot = lambda7;
lambda1dot = -0.5*(-lambda2*w3 + lambda3*w2 - lambda4*w1);
lambda2dot = -0.5*(lambda1*w3 - lambda4*w2 - lambda3*w1);
lambda3dot = -0.5*(-lambda4*w3 - lambda1*w2 + lambda2*w1);
lambda4dot = -0.5*(lambda3*w3 + lambda2*w2 + lambda1*w1);
lambda5dot = -0.5*(lambda1*q4 + lambda2*q3 - lambda3*q2 - lambda4*q1);
lambda6dot = -0.5*(lambda2*q4 - lambda1*q3 - lambda4*q2 + lambda3*q1);
lambda7dot = -0.5*(lambda3*q4 + lambda1*q2 - lambda4*q3 - lambda2*q1);


states = [q1dot q2dot q3dot q4dot w1dot w2dot w3dot lambda1dot lambda2dot lambda3dot lambda4dot lambda5dot lambda6dot lambda7dot]';

