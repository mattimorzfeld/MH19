function [x,resnorm,residual,exitflag,output,lambda,jacobian]=myMinLS2(xo,z,dT,dt,H,R,mu,Lb)
func=@(x)funcF2(x,z,dT,dt,H,R,mu,Lb);
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
    'Display','off');%'iter-detailed');
n = size(H,2);
ub = inf*ones(n,1);
lb = [-inf;-inf;0];
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(func,xo,lb,ub,options);


