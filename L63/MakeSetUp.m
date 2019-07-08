%%
clear
close all
clc
%% ------------------------------------------------------------------------


%% ------------------------------------------------------------------------
dt = .05;
nAssims = 1200;
SetUp = 4; 
noisy = 0;
sqrtQ = noisy*sqrt(dt)*eye(3);

%% Background covariance
xo = [-2.7295;-0.6538;24.2952];
[yC,~]=RunModel(xo,dt,10000,sqrtQ);
B = cov(yC');
sqrtB = sqrtm(B);
       
for Gap = [2 4 6 8 10 12 14]
    dT = Gap*dt;
    T = nAssims*Gap*dt;
    
    %% Simulation of the truth
    [y,t]=RunModel(xo,dt,T,sqrtQ);

    %% Observation network
    if SetUp == 1
        H = [1 0 0;0 1 0; 0 0 1];
        nObs = size(H,1);
        R = 0.1*ones(nObs,nAssims);
    elseif SetUp == 2
        H = [1 0 0;0 1 0; 0 0 1];
        nObs = size(H,1);
        R = ones(nObs,nAssims);
    elseif SetUp == 3
        H = [1 0 0;0 1 0; 0 0 1];
        nObs = size(H,1);
        R = 2*ones(nObs,nAssims);
    elseif SetUp == 4
        H = [1 0 0; 0 0 1];
        nObs = size(H,1);
        R = 0.1*ones(nObs,nAssims);
    end
    
    %% get obs
    z = zeros(nObs,nAssims);
    tObs = zeros(1,nAssims);
    for kk=1:nAssims
        z(:,kk) = H*y(:,kk*Gap+1);
        z(:,kk) = z(:,kk)+sqrt(R(:,kk)).*randn(length(z(:,kk)),1);
        tObs(kk) = t(kk*Gap+1);
    end
    
    %% Save
    FileName = strcat('./SetUps/SetUp_',num2str(SetUp),'_Gap_',num2str(Gap),'.mat');
    save(FileName,'H','R','dt','dT','Gap','nAssims','nObs','B','sqrtB','sqrtQ','y','t','z','tObs','xo');
end
%% ------------------------------------------------------------------------