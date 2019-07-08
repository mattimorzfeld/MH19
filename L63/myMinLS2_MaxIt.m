function [x,resnorm,residual,exitflag,output,lambda,jacobian]=myMinLS2_MaxIt(xo,z,dT,dt,H,R,mu,Lb,MaxIt)
func=@(x)funcF2(x,z,dT,dt,H,R,mu,Lb);
options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective',...
    'MaxIteration',MaxIt, ...
    'Display','off');%'iter-detailed');
n = size(H,2);
ub = inf*ones(n,1);
lb = [-inf;-inf;0];
[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(func,xo,lb,ub,options);


