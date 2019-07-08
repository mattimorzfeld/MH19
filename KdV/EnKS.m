function [xEnKF,tracePEnKF,X1a]=EnKS(X,Z,locrad,infl,H,R,Nx,Nobs,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt)

Ne = size(X,2);
xEnKF = zeros(Nx,nAssims);
tracePEnKF = zeros(1,nAssims); 

dsm = locrad; % Localization tuning
build_localization

X0 = X;
X0a = X0;
X1 = X0;
X1a = X0;
for kk=1:nAssims
    fprintf('EnKS assim %g / %g\r',kk,nAssims)
    % forecast
    for jj=1:Ne
        Cv = run_model_RK3(X0(:,jj), Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        X1(:,jj) = tmp(:,end);
    end

    % EnKS
    X0m = mean(X0,2);
    X0pert = X0 - X0m*ones(1,Ne);
    X0 = X0m*ones(1,Ne)+sqrt(infl)*X0pert;
    X1m = mean(X1,2);
    X1pert = X1 - X1m*ones(1,Ne);
    X1 = X1m*ones(1,Ne)+sqrt(infl)*X1pert;
    
    P01 = CL.*(X0pert*X1pert'/(Ne-1)); % prior covariance from time 1 to time 0
    P11 = CL.*(X1pert*X1pert'/(Ne-1)); % prior covariance from time 1 to time 1
    K = P01*H'/(H*P11*H'+diag(R(:,kk))); 
    
    for jj = 1:Ne
        X0a(:,jj) = X0(:,jj)+K*(Z(:,kk)+sqrtm(diag(R(:,kk)))*randn(Nobs,1)-H*X1(:,jj));
        
        % push analysis ensemble at time 0 to time 1
        Cv = run_model_RK3(X0a(:,jj), Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        X1a(:,jj) = tmp(:,end);
    end
    % handle mean
    X0ma = X0m+K*(Z(:,kk)-H*X1m);
    Cv = run_model_RK3(X0ma, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
    tmp = ifft( Cv.','Symmetric' );
    X1ma = tmp(:,end);
            
    % save results
    xEnKF(:,kk) = X1ma;
    tracePEnKF(kk) = trace(cov(X1a'));
    % reset ensemble (make time 1 time 0)
    X0 = X1a;
end