function [Z,Xt,R] = SynthExp(xt,r,H,Nobs,nAssims,Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt)

Xt = zeros(Nx,nAssims);
Z = zeros(Nobs,nAssims);
R = zeros(Nobs,nAssims);
for kk=1:nAssims
    Cv = run_model_RK3(xt, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
    tmp = ifft( Cv.','Symmetric' );
    xtT = tmp(:,end);
    R(:,kk) = r*abs(H*xtT);
    z = H*xtT+sqrt(R(:,kk)).*randn(Nobs,1);
    xt = xtT;
    
    % save results
    Z(:,kk) = z;
    Xt(:,kk) = xt;
end