function [Xa,xAll,traceP] = PF_Pot(Ne,a,r,z,xo,sqrtPo,dt,dT,sqrtQ,H,R)

nAssims = size(z,2);
nSteps = floor(dT/dt*nAssims+1);
Gap = floor(dT/dt);

%% initial ensemble
X = zeros(3,Ne);
xpf = cell(Ne,1);
for ll=1:Ne
    X(:,ll) = xo+sqrtPo*randn(3,1);
    xpf{ll} = X(:,ll)';
end

%% cycle
xAll = zeros(3,nSteps);
Xa = zeros(3,nAssims);
traceP = zeros(1,nAssims);
for kk=1:nAssims
    fprintf('Poterjoys PF Assim %g / %g\r',kk,nAssims)
    RunningMean = zeros(3,Gap+1);
    for ll=1:Ne
        tmp = RunModel(X(:,ll),dt,dT,sqrtQ);
        RunningMean = RunningMean+tmp;
        xpf{ll} = tmp(:,end)';
    end
    RunningMean = RunningMean/Ne;
    
    
    [~,xpf,eflag] = pf_update_MWR16_beta(xpf,3,Ne,H,z(:,kk)',r,a,1:size(H,1),R(:,kk),0);
    if eflag == 0
        fprintf('Error in localized PF.\n')
    end

    %% save results
    for oo=1:Ne
        X(:,oo) = xpf{oo}';
    end  

    Xa(:,kk) = mean(X,2);
    traceP(kk) = trace(cov(X'));
    xAll(:,(kk-1)*Gap+1:kk*Gap+1)=RunningMean;
    xAll(:,kk*Gap+1)=Xa(:,kk);
end