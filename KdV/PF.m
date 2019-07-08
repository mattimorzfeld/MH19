function [xPF,tracePPF,Xa,rho]=PF(X,Z,locrad,infl,H,R,Nx,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt)

Ne = size(X,2);
xPF = zeros(Nx,nAssims);
tracePPF = zeros(1,nAssims); 
rho = zeros(nAssims,1);

dsm = locrad; % Localization tuning
build_localization

Xa = zeros(Nx,Ne);
for kk=1:nAssims
    fprintf('PF assim %g / %g',kk,nAssims)
    % forecast
    w = zeros(Ne,1);
    for jj=1:Ne
        Cv = run_model_RK3(X(:,jj), Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        X(:,jj) = tmp(:,end);
        w(jj) = norm(sqrt(2*R(:,kk)).\(Z(:,kk)-H*X(:,jj)))^2;
    end
    w = normalizeweights(w);
    rho(kk) = mean(w.^2)/mean(w)^2;
    fprintf('    Effective sample size: %g\r',Ne/rho(kk))
    
    ma = X*w;
    Xa =  resampling(w,X,Ne,Nx);
    Pa = infl*CL.*cov(Xa');
    sqrtPa = chol(Pa)';
    X = MakeInitialEnsemble(ma,sqrtPa,Ne,Nx);
    
    
    %% save results
    xPF(:,kk) = ma;
    tracePPF(kk) = trace(Pa);
end