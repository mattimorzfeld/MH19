function [y,t]=RunModel(xo,dt,T,sqrtQ)

% [t,y] = ode45(@nm,[0:dt:T],xo);
% y = y';

t = 0:dt:T;
Steps = length(t);
y = zeros(3,Steps);
y(:,1) = xo;

for kk = 2:Steps    
     k1 = L63(t(kk),y(:,kk-1));
     k2 = L63(t(kk)+dt/2,y(:,kk-1)+dt/2*k1);
     k3 = L63(t(kk)+dt/2,y(:,kk-1)+dt/2*k2);
     k4 = L63(t(kk)+dt,y(:,kk-1)+dt*k3);
     RK4 = y(:,kk-1)+dt/6*(k1+2*k2+2*k3+k4);
     y(:,kk) = RK4+sqrt(dt)*sqrtQ*randn(3,1);
end

