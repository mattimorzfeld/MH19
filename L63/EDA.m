function [Xa,xAll,traceP] = EDA(Ne,z,mu,sqrtB,infl,dt,dT,H,R)

nAssims = size(z,2);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);

%% initial ensemble
X = zeros(3,Ne);
for ll=1:Ne
    X(:,ll) = mu+sqrtB*randn(3,1);
end

%% cycle
xAll = zeros(3,nSteps);
Xa = zeros(3,nAssims);
traceP = zeros(1,nAssims);

for kk=1:nAssims
    fprintf('EDA Assim %g / %g\r',kk,nAssims)
    % forecast
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        tmp = RunModel(X(:,ll),dt,dT,0*eye(3));
        RunningMean = RunningMean+tmp;
    end
    RunningMean = RunningMean/Ne;
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    
    % data assimilation: ensemble 4d-Var
    Xo = zeros(3,Ne);
    Xt = zeros(3,Ne);
    [x0star,~,~,~,~,~,~]=myMinLS2(mu,z(:,kk), ...
                                dT,dt,H,R(:,kk),mu,sqrtB); 
    Xo(:,1)= x0star; 
    [tmp,~] = RunModel(Xo(:,1),dt,dT,0*eye(3));
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=tmp;
    Xa(:,kk) = tmp(:,end);
    Xt(:,1) = tmp(:,end);
    % ensemble
    for ll = 2:Ne
        [x0star,~,~,~,~,~,~]=myMinLS2(mu,z(:,kk)+sqrt(R(:,kk)).*randn(size(H,1),1), ...
                                dT,dt,H,R(:,kk),mu+sqrtB*randn(3,1),sqrtB); 
        Xo(:,ll)= x0star;
        [tmp,~] = RunModel(x0star,dt,dT,0*eye(3));
        Xt(:,ll) = tmp(:,end);
    end
   
    traceP(kk) = trace(cov(Xt'));  

    % data assimilation: update background
    mu = Xt(:,1);
    sqrtB = infl*sqrtm(cov(Xt'));
end