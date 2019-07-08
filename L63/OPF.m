function [Xa,xAll,rho,traceP] = OPF(Ne,z,xo,sqrtPo,dt,dT,sqrtQ,H,R)

nAssims = size(z,2);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);
Q = dt*sqrtQ*sqrtQ;

%% initial ensemble
X = zeros(3,Ne);
rho  = zeros(nAssims,1);
for ll=1:Ne
    X(:,ll) = constGauss(xo,sqrtPo);
end

%% cycle
xAll = zeros(3,nSteps);
Xa = zeros(3,nAssims);
traceP = zeros(1,nAssims);
for kk=1:nAssims
    fprintf('OPF Assim %g / %g\r',kk,nAssims)
    w = zeros(Ne,1);
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        tmp = RunModel(X(:,ll),dt,dT,sqrtQ);
        RunningMean = RunningMean+tmp;
        tmp = RunModel(tmp(:,end-1),dt,dt,0*sqrtQ);
        fxnm1 = real(tmp(:,end));
        
        Rtmp = diag(R(:,kk));
        K = (Q*H')/(Rtmp+H*Q*H');
        sig = (eye(3)-K*H)*Q;
%         sig
%         Rtmp
%         Q
%         z(:,kk)
%         Rtmp\z(:,kk)
%         Q\fxnm1+Rtmp\z(:,kk)
% z(:,kk)
% Rtmp
%  Q\fxnm1+H'*(Rtmp\z(:,kk))
        tmp = sig*(Q\fxnm1 + H'*(Rtmp\z(:,kk)));
        X(:,ll) = tmp+sqrtm(sig)*randn(3,1);
        w(ll) = .5*norm(sqrtm(H*Q*H'+Rtmp)\(z(:,kk)-H*fxnm1))^2;
          
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