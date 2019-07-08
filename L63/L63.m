function dydt = nm(t,x)
sigma = 10; beta=8/3; rho = 28;
dydt = [sigma* (x(2)-x(1)); ...
        x(1)*(rho-x(3))-x(2);
        x(1)*x(2)-beta*x(3)];
% al = .3;
% be =-1;
% gam = 1;
% del = .4;
% om = 1.2;
% % x(3) = x(3)+.18;
% 
% dydt = [x(2); 
%                 -al*x(2)-be*x(1)-gam*x(1)^3+del*cos(om*t);
%                 2*sqrt(x(3)+0.25)*x(2)-x(2)-1e-1*(x(3)-0.006)
%                 ];