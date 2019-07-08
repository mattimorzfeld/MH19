function [xEnKF,tracePEnKF,Xa]=EnKF(X,Z,locrad,infl,H,R,Nx,Nobs,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt)

Ne = size(X,2);
xEnKF = zeros(Nx,nAssims);
tracePEnKF = zeros(1,nAssims); 

dsm = locrad; % Localization tuning
build_localization

Xa = zeros(Nx,Ne);
for kk=1:nAssims
    fprintf('EnKF assim %g / %g\r',kk,nAssims)
    % forecast
    for jj=1:Ne
        Cv = run_model_RK3(X(:,jj), Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        X(:,jj) = tmp(:,end);
    end
    mf = mean(X,2);
    Pf = infl*(CL.*cov(X'));

    K = (Pf*H')/(H*Pf*H'+diag(R(:,kk))*eye(Nobs));
    Xpert = X - mf*ones(1,Ne);
    X = mf*ones(1,Ne)+sqrt(infl)*Xpert;
    
    for jj = 1:Ne
        Xa(:,jj) = X(:,jj)+K*(Z(:,kk)+sqrtm(diag(R(:,kk)))*randn(Nobs,1)-H*X(:,jj));
    end
    ma = mf+K*(Z(:,kk)-H*mf);
    X = Xa;
        
    %% save results
    xEnKF(:,kk) = ma;
    tracePEnKF(kk) = trace(cov(Xa'));
end