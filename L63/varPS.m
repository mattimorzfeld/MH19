function [Xa,xAll,rho,traceP] = varPS(Ne,z,mu,sqrtB,infl,dt,dT,H,R)

nAssims = size(z,2);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);

%% initial ensemble
X = zeros(3,Ne);
rho  = zeros(nAssims,1);
for ll=1:Ne
    X(:,ll) = mu+sqrtB*randn(3,1);
end

%% cycle
xAll = zeros(3,nSteps);
Xa = zeros(3,nAssims);
traceP = zeros(1,nAssims);

for kk=1:nAssims
    fprintf('varPS Assim %g / %g\r',kk,nAssims)
    % forecast
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        tmp = RunModel(X(:,ll),dt,dT,0*eye(3));
        RunningMean = RunningMean+tmp;
    end
    RunningMean = RunningMean/Ne;
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    
    % data assimilation: 4d-Var
    [x0star,~,~,~,~,~,J]=myMinLS2(mu,z(:,kk),dT,dt,H,R(:,kk),mu,sqrtB); 
    Co = (2*(J'*J))\eye(3);
    sqrtCo = sqrtm(Co);   
    
    w = zeros(1,Ne);
    Xt = zeros(3,Ne);    
        
    [tmp,~] = RunModel(x0star,dt,dT,0*eye(3));
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=tmp;
    Xa(:,kk) = tmp(:,end);
    
    % data assimilation: ensemble generation
    Xt(:,1) = tmp(:,end);
    f =  funcF2(x0star,z(:,kk),dT,dt,H,R(:,kk),mu,sqrtB);
    F = f'*f;
    w(1) = F;
    for ll=2:Ne
        xo = x0star+sqrtCo*randn(3,1);
        tmp = RunModel(xo,dt,dT,0*eye(3));
        Xt(:,ll) = tmp(:,end);
        f =  funcF2(xo,z(:,kk),dT,dt,H,R(:,kk),mu,sqrtB);
        F = f'*f;
        Fo = .5*norm(sqrtCo\(xo-x0star))^2;
        w(ll) = F-Fo;
    end
    w = normalizeweights(w);
    Xt = resampling(w,Xt,Ne,3);
    rho(kk) = mean(w.^2)/mean(w)^2;
    traceP(kk) = trace(cov(Xt'));  

    % data assimilation: update background
    mu = Xt(:,1);
    sqrtB = infl*sqrtm(cov(Xt'));
    for jj=1:3
        if sqrtB(jj,jj)<sqrt(.001*abs(Xt(jj,1)))
            sqrtB(jj,jj) = sqrt(.001*abs(Xt(jj,1)));
        end
    end
end