function [Xa,xAll,rho,traceP,Xrs] = PF_GR_saveNonGaussFeatures(Ne,infl,z,xo,sqrtPo,dt,dT,sqrtQ,H,R)

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

FileName = strcat('PFGR_SaveNonGauss_dT_',num2str(dT),'_infl_',num2str(infl),'_Ne_',num2str(Ne),'.mat');

ForecastKurt = zeros(nAssims,3);
ForecastSkew = zeros(nAssims,3);

AnAt0Kurt = zeros(nAssims,3);
AnAt0Skew = zeros(nAssims,3);

AnAt1Kurt = zeros(nAssims,3);
AnAt1Skew = zeros(nAssims,3);

for kk=1:nAssims
    fprintf('PF w GR Assim %g / %g\r',kk,nAssims)
    Xo = X; % store ensemble prior to analysis time
    w = zeros(Ne,1);
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        tmp = RunModel(X(:,ll),dt,dT,sqrtQ);
        RunningMean = RunningMean+tmp;
        X(:,ll) = real(tmp(:,end));
        w(ll) = .5*(z(:,kk)-H*X(:,ll))'*(R(:,kk).\(z(:,kk)-H*X(:,ll)));
    end
    %% Skewness and kurtosis of forecast ensemble
    ForecastKurt(kk,:) = kurtosis(X',1);
    ForecastSkew(kk,:) = skewness(X');
    
    RunningMean = RunningMean/Ne;
    w = normalizeweights(w);
    rho(kk) = mean(w.^2)/mean(w)^2;
    Xrs = resampling(w,X,Ne,3);
    
    %% Skewness and kurtosis of analysis ensemble at obs time
    AnAt1Kurt(kk,:) = kurtosis(Xrs',1);
    AnAt1Skew(kk,:) = skewness(Xrs');
    
    %% Skewness and kurtosis of analysis ensemble at time 0
    Xors = resampling(w,Xo,Ne,3);
    AnAt0Kurt(kk,:) = kurtosis(Xors',1);
    AnAt0Skew(kk,:) = skewness(Xors');
    
    Xa(:,kk) = X*w;
    traceP(kk) = trace(cov(Xrs'));
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    xAll(:,kk*Gap+1)=Xa(:,kk);
    
    %% draw new ensemble
    sqrtC = infl*sqrtm(cov(Xrs'));
    XApp = zeros(3,Ne);
    for ll=2:Ne
        XApp(:,ll) = Xa(:,kk)+sqrtC*randn(3,1);
    end
    X = XApp;
    
end
save(FileName,'AnAt0Kurt','AnAt0Skew','AnAt1Kurt','AnAt1Skew','ForecastKurt','ForecastSkew')
