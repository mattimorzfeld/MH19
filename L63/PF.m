function [Xa,xAll,rho,traceP] = PF(Ne,z,xo,sqrtPo,dt,dT,sqrtQ,H,R)

nAssims = size(z,2);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);

%% initial ensemble
X = zeros(3,Ne);
rho  = zeros(nAssims,1);
for ll=1:Ne
    X(:,ll) = xo+sqrtPo*randn(3,1);
end

%% cycle
xAll = zeros(3,nSteps);
Xa = zeros(3,nAssims);
traceP = zeros(1,nAssims);
for kk=1:nAssims
    fprintf('PF Assim %g / %g\r',kk,nAssims)
    w = zeros(Ne,1);
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        tmp = RunModel(X(:,ll),dt,dT,sqrtQ);
        RunningMean = RunningMean+tmp;
        X(:,ll) = real(tmp(:,end));
        w(ll) = .5*(z(:,kk)-H*X(:,ll))'*(R(:,kk).\(z(:,kk)-H*X(:,ll)));
    end
    RunningMean = RunningMean/Ne;
    w = normalizeweights(w);
    rho(kk) = mean(w.^2)/mean(w)^2;
    Xrs = resampling(w,X,Ne,3);

    Xa(:,kk) = X*w;
    traceP(kk) = trace(cov(Xrs'));
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    xAll(:,kk*Gap+1)=Xa(:,kk);
    
    X = Xrs;
end