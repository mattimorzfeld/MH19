function [xPF,tracePPF,Xa]=lPF(X,Z,locrad,infl,H,R,Nx,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt,ob_step,HowManyObs)

Ne = size(X,2);
xPF = zeros(Nx,nAssims);
tracePPF = zeros(1,nAssims); 
rho = zeros(nAssims,1);

nObs = size(H,1);

dsm = locrad; % Localization tuning
build_localization

Xa = zeros(Nx,Ne);
for kk=1:nAssims
    fprintf('l-PF assim %g / %g\r',kk,nAssims)
    % forecast
    for jj=1:Ne
        Cv = run_model_RK3(X(:,jj), Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        X(:,jj) = tmp(:,end);
    end

    for oo=1:Nx
        w = zeros(Ne,1);
        WhichObs = ceil(oo/ob_step)+(-HowManyObs:HowManyObs);
        inds = find(WhichObs>0);
        WhichObs = WhichObs(inds);
        inds = find(WhichObs<=nObs);
        WhichObs = WhichObs(inds);
        for jj=1:Ne
            w(jj) = norm(sqrt(2*R(WhichObs,kk)).\(Z(WhichObs,kk)-H(WhichObs,:)*X(:,jj)))^2;
        end
        w = normalizeweights(w);
        Xa(oo,:) = resamplingScalar(w,X(oo,:),Ne);
    end
    ma = mean(Xa,2);
    Pa = infl*(CL.*cov(Xa'));
    [U,L] = eig(Pa);
    sqrtPa = U*sqrt(abs(L));
    X = MakeInitialEnsemble(ma,sqrtPa,Ne,Nx);
    
    %% save results
    xPF(:,kk) = ma;
    tracePPF(kk) = trace(Pa);%sum(sum(abs(L)));
end
