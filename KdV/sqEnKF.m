function [xEnKF,tracePEnKF,Xa]=sqEnKF(X,Z,locrad,infl,H,R,Nx,Nobs,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt)

Ne = size(X,2);
xEnKF = zeros(Nx,nAssims);
tracePEnKF = zeros(1,nAssims); 

dsm = locrad; % Localization tuning
build_localization

Xa = zeros(Nx,Ne);
for kk=1:nAssims
    fprintf('sqEnKF assim %g / %g\r',kk,nAssims)
    % forecast
    for jj=1:Ne
        Cv = run_model_RK3(X(:,jj), Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        X(:,jj) = tmp(:,end);
    end
    mf = mean(X,2);
    Pf = infl*(CL.*cov(X'));

    K = (Pf*H')/(H*Pf*H'+diag(R(:,kk))*eye(Nobs));

    ma = mf+K*(Z(:,kk)-H*mf);
    Xpert = X - mf*ones(1,Ne);
    X = mf*ones(1,Ne)+sqrt(infl)*Xpert;
    
    z = (sqrt(infl)/sqrt(Ne-1))*Xpert;

    tmp = z'*H'*(R(:,kk).\(H*z));
    tmp = .5*(tmp+tmp');
    [E,OM] = eig(tmp);
    T = E*(sqrtm(eye(Ne)+OM)\E');
    Za = z*T;
    X = repmat(ma,1,Ne)+sqrt(Ne-1)*Za;
        
    %% save results
    xEnKF(:,kk) = ma;
    tracePEnKF(kk) = trace(cov(X'));
end