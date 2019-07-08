function [XamSave,xAll,traceP,X] = mySqEnKF(Ne,infl,z,dt,dT,sqrtQ,H,R)
nAssims = size(z,2);
nObs = size(z,1);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);

%% pre-allocate
traceP = zeros(nAssims,1);
D = zeros(nObs,Ne);
xAll = zeros(3,nSteps);
XamSave = zeros(3,nAssims);

%% initial ensemble
xo = [0.4424    0.3260    0.0308]';
[yC,~]=RunModel(xo,dt,10000,sqrtQ);
X = yC(:,randi(length(yC),Ne,1));

for kk=1:nAssims
    fprintf('SqEnKF assim %g / %g\r',kk,nAssims)
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        trajectory = RunModel(X(:,ll),dt,dT,sqrtQ);
        X(:,ll)=trajectory(:,end);
        RunningMean=RunningMean+trajectory;
        D(:,ll) = z(:,kk)+sqrt(R(:,kk)).*randn(nObs,1);
    end 
    % Intermediate time steps
    RunningMean = RunningMean/Ne;
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;

    P = infl*cov(X');
    
    Xm = mean(X,2);
    Xpert = X - Xm*ones(1,Ne);
    
    K = P*H'*((H*P*H'+diag(R(:,kk)))\eye(nObs));
    Xam = Xm + K*(z(:,kk)-H*Xm);
    xAll(:,kk*Gap+1)= Xam;
   
    Z = (sqrt(infl)/sqrt(Ne-1))*Xpert;

    tmp = Z'*H'*(R(:,kk).\(H*Z));
    tmp = .5*(tmp+tmp');
    [E,OM] = eig(tmp);
    T = E*(sqrtm(eye(Ne)+OM)\E');
    Za = Z*T;
    X = repmat(Xam,1,Ne)+sqrt(Ne-1)*Za;
    
    XamSave(:,kk) = Xam;
    traceP(kk) = trace(cov(X'));
end