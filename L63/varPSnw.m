function [Xa,xAll,traceP] = varPSnw(Ne,z,mu,sqrtB,infl,dt,dT,H,R)

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
    fprintf('varPS no weights Assim %g / %g\r',kk,nAssims)
    % forecast
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        tmp = RunModel(X(:,ll),dt,dT,0*eye(3));
        RunningMean = RunningMean+tmp;
    end
    RunningMean = RunningMean/Ne;
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    
    % data assimilation: 4D-Var
    [x0star,~,~,~,~,~,J]=myMinLS2(mu,z(:,kk),dT,dt,H,R(:,kk),mu,sqrtB); 
    Co = (2*(J'*J))\eye(3);
    sqrtCo = sqrtm(Co);   
    
    Xt = zeros(3,Ne);    
    tmp = RunModel(x0star,dt,dT,0*eye(3));
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=tmp;
    Xa(:,kk) = tmp(:,end);
    Xt(:,1) = tmp(:,end);
    
    % data assimilation: ensemble generation
    for ll=2:Ne
        xo = x0star+sqrtCo*randn(3,1);
        tmp = RunModel(xo,dt,dT,0*eye(3));
        Xt(:,ll) = tmp(:,end);
    end
    traceP(kk) = trace(cov(Xt'));  

    % data assimilation: background update
    mu = Xt(:,1);
    sqrtB = infl*sqrtm(cov(Xt'));
%     for jj=1:3
%         if sqrtB(jj,jj)<sqrt(.001*abs(Xt(jj,1)))
%             sqrtB(jj,jj) = sqrt(.001*abs(Xt(jj,1)));
%         end
%     end
end