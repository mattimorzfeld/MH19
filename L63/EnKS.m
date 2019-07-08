function [Xam,xAll,traceP] = EnKS(Ne,infl,z,dt,dT,sqrtQ,H,R)

nAssims = size(z,2);
nObs = size(z,1);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);

%% pre-allocate
traceP = zeros(nAssims,1);
D = zeros(nObs,Ne);
xAll = zeros(3,nSteps);
X0am = zeros(3,nAssims);
Xam = zeros(3,nAssims);

%% initial ensemble
xo = [0.4424    0.3260    0.0308]';
[yC,~]=RunModel(xo,dt,10000,sqrtQ);
X0 = yC(:,randi(length(yC),Ne,1));
X1 = X0; % just initialize X1 here ...
%% cycle
for kk=1:nAssims
    fprintf('EnKS Assim %g / %g\r',kk,nAssims)
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        trajectory = RunModel(X0(:,ll),dt,dT,sqrtQ);
        X1(:,ll)=trajectory(:,end);
        RunningMean=RunningMean+trajectory;
        D(:,ll) = z(:,kk)+sqrt(R(:,kk)).*randn(nObs,1);
    end 
    % Intermediate time steps
    RunningMean = RunningMean/Ne;
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    
    % EnKS
    X0m = mean(X0,2);
    X0pert = X0 - X0m*ones(1,Ne);
    X0 = X0m*ones(1,Ne)+sqrt(infl)*X0pert;
    X1m = mean(X1,2);
    X1pert = X1 - X1m*ones(1,Ne);
    X1 = X1m*ones(1,Ne)+sqrt(infl)*X1pert;
    
    P01 = X0pert*X1pert'/(Ne-1); % prior covariance from time 1 to time 0
    P11 = X1pert*X1pert'/(Ne-1); % prior covariance from time 1 to time 1
    
    K = P01*H'/(H*P11*H'+diag(R(:,kk))); % Kalman gain for observation at time 1 but update at time 0
   
    X0a = X0+K*(D-H*X1); % perturbed obs
    X0am(:,kk) = X0m + K*(z(:,kk)-H*X1m); % analysis
    
    % push posterior ensemble at time 0 to time 1, and re-label to time 0 to restart loop
    X0 = X0a;
    for ll=1:Ne
        trajectory = RunModel(X0(:,ll),dt,dT,sqrtQ);
        X0(:,ll) = trajectory(:,end); 
    end
    
    % get mean
    trajectory = RunModel(X0am(:,kk),dt,dT,sqrtQ);
    Xam(:,kk) = trajectory(:,end); 
    
    % I think you want to save your stats at time 1, but I'm not sure ... 
    traceP(kk) = trace(cov(X0')); 
    xAll(:,kk*Gap+1) = mean(X0,2); 
end