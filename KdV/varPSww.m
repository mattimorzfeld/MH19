function [xvarPS,tracePvarPS,Xa_varPS]=varPSww(Ne,mu,Lb,Z,locrad,infl,H,R,Nx,nAssims, ...
                                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt, OuterLoops)

xvarPS = zeros(Nx,nAssims);
tracePvarPS = zeros(1,nAssims); 

dsm = locrad; % Localization tuning
build_localization
Xa_varPS = zeros(Nx,Ne);
for kk=1:nAssims
    fprintf('varPS assim %g / %g\r',kk,nAssims)
    [mx4DVar,~,~,~,~,~,J] = ...
                myMinLS2(zeros(Nx,1),Z(:,kk),H,R(:,kk),mu,Lb, ...
                            Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt, OuterLoops);
    x4DVar = mu+Lb*mx4DVar;   
    
    B4DVar = (2*(J'*J))\eye(Nx);
    B4DVar = CL.*B4DVar;
    LB4DVar = chol(B4DVar)';
    Ens4DVar = zeros(Nx,Ne);
    Ens4DVar(:,1) = x4DVar;
    for jj=2:Ne
        xi = Lb*LB4DVar*randn(Nx,1);
        go = 1;
        while go == 1
            xi = Lb*(LB4DVar*randn(Nx,1));
            sp = x4DVar + xi;
            sm = x4DVar - xi;
            if min(sp)>-1
%                 fprintf('Found 4DVar Ens %g/ %g\n',jj,Ne)
                Ens4DVar(:,jj) = sp;
                go = 0;
            elseif min(sm)>-1
%                 fprintf('Found 4DVar Ens %g/ %g\n',jj,Ne)
                Ens4DVar(:,jj) = sm;
                go = 0;
            end
        end
    end

    for jj = 1:Ne
        Cv = run_model_RK3(Ens4DVar(:,jj), Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        Xa_varPS(:,jj) = tmp(:,end);
    end
    % update background
    mu = Xa_varPS(:,1);
    B = infl*CL.*cov(Xa_varPS');
    Lb = chol(B)';
    % save
    tracePvarPS(kk) = trace(B);
    xvarPS(:,kk) = mu;

    
end