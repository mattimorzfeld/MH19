function [f,J] = funcF2(x,z,H,R,mu,Lb,Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt)
k = size(H,1); % number of obs

X = mu+Lb*x;

Cv = run_model_RK3(X, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
trajectory = ifft( Cv.','Symmetric' );
X = trajectory(:,end);
f = zeros(k,1);
for kk=1:k
    f(kk) = (H(kk,:)*X-z(kk))./sqrt(2*R(kk));
end
f = real([sqrt(.5)*x; f]);

if nargout >1
    % build matrix TLM from model
    M = zeros(Nx,Nx);
    parfor j = 1:Nx
        vect = zeros(Nx,1);
        vect(j) = 1;
        CT = run_TLM(vect, trajectory, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        mw = ifft(CT(Nt,:),'symmetric');
        M(:,j) = mw;
    end
    J = [sqrt(.5)*eye(Nx); sqrt(2*R).\(H*M)];
end