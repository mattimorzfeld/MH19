function [xAll,traceP]=EDA(Ne,mu,Lb,Z,locrad,infl,H,R,Nx,nAssims, ...
                                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt, OuterLoops)

xAll = zeros(Nx,nAssims);
traceP = zeros(1,nAssims); 

dsm = locrad; % Localization tuning
build_localization
X = zeros(Nx,Ne);
for kk=1:nAssims
    fprintf('EDA assim %g / %g\n',kk,nAssims)
    [mx4DVar,~,~,~,~,~,~] = ...
        myMinLS2(zeros(Nx,1),Z(:,kk),H,R(:,kk),mu,Lb, ...
        Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt,OuterLoops);
    x4DVar = mu+Lb*mx4DVar;
    
    Cv = run_model_RK3(x4DVar, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
    tmp = ifft( Cv.','Symmetric' );
    x4DVarT = tmp(:,end);
    X(:,1) = x4DVarT;
    
    % ensemble
    for jj=2:Ne
        [mx4DVar,~,~,~,~,~,~] = ...
            myMinLS2(zeros(Nx,1),Z(:,kk)+sqrt(R(:,kk)).*randn(size(H,1)),H,R(:,kk),mu+Lb*randn(Nx,1),Lb, ...
            Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt, OuterLoops);
        x4DVar = mu+Lb*mx4DVar;
        Cv = run_model_RK3(x4DVar, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        X(:,jj) = tmp(:,end);
    end
    Pa = infl*CL.*cov(X');
    [U,L] = eig(Pa);
    sqrtPa = U*sqrt(abs(L));

    % save
    traceP(kk) = sum(sum(abs(L)));
    xAll(:,kk) = X(:,1);
    % update background
    mu = X(:,1);
    B = U*sqrt(abs(L));
end