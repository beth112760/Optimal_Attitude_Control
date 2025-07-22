function phid = rbf1(c,t,ti)


r = t - ti;
R = abs(r);

h = 1/c*R;

%Multiquadric
phid = (h./c./sqrt(h.^2+1) + 5*R.^4).*sign(r);


%Inverse Multiquadric
%phid = (-c.*R./(c.^2+R.^2).^(3/2)+ 5*R.^4).*sign(r);
%Gaussian
%phid = (-(2*exp(-h)).*h + 5*R.^4).*sign(r);
%Inverse quadratic
%phid = (-2*h./((h.^2+1).^(3/2).*sqrt(h.^2+1))+ 5*R.^4).*sign(r);

end