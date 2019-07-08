function [Xam,xAll,traceP] = EnKF(Ne,infl,z,dt,dT,sqrtQ,H,R)

nAssims = size(z,2);
nObs = size(z,1);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);

%% pre-allocate
traceP = zeros(nAssims,1);
D = zeros(nObs,Ne);
xAll = zeros(3,nSteps);
Xam = zeros(3,nAssims);

%% initial ensemble
xo = [0.4424    0.3260    0.0308]';
[yC,~]=RunModel(xo,dt,10000,sqrtQ);
X = yC(:,randi(length(yC),Ne,1));

%% cycle
for kk=1:nAssims
    fprintf('EnKF Assim %g / %g\r',kk,nAssims)
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
    
    % EnKF
    P = infl*cov(X');
    
    Xm = mean(X,2);
    Xpert = X - Xm*ones(1,Ne);
    X = Xm*ones(1,Ne)+sqrt(infl)*Xpert;
    K = P*H'/(H*P*H'+diag(R(:,kk)));
   
    Xa = X+K*(D-H*X);
    Xam(:,kk) = Xm + K*(z(:,kk)-H*Xm);
    traceP(kk) = trace(cov(Xa'));
        
    X = Xa;
    xAll(:,kk*Gap+1) = Xam(:,kk) ;
end