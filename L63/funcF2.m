function f = funcF2(x,z,dT,dt,H,R,mu,Lb)
trajectory = RunModel(x,dt,dT,0*eye(3));
X = trajectory(:,end);
k = size(H,1); % number of obs
f = zeros(k,1);
for kk=1:k
    f(kk) = (H(kk,:)*X-z(kk))/sqrt(2*R(kk));
end
f = real([sqrt(.5)*(Lb\(x-mu)); f]);


