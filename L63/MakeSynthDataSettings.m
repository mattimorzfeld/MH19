%%
clear
close all
clc
%% ------------------------------------------------------------------------


% %% Set up
% %% ------------------------------------------------------------------------
% xo = [0.4424    0.3260    0.0308]';
% 
% dt = .05;
% Gap = 13;
% dT = Gap*dt;
% nAssims = 1500;
% T = nAssims*Gap*dt;
% sqrtQ = 1*sqrt(dt)*diag([1 1 1]);
% 
% %% Background covariance
% xo = [-2.7295;-0.6538;24.2952];
% [yC,~]=RunModel(xo,dt,10000,sqrtQ);
% B = cov(yC');
% sqrtB = sqrtm(B);
% 
% [y,t]=RunModel(xo,dt,T,sqrtQ);
% 
% H = [1 0 0;0 1 0; 0 0 1];
% MinNoise = [.1;.1;.01];
% 
% nObs = size(H,1);
% 
% R = zeros(nObs,nAssims);
% z = zeros(nObs,nAssims);
% tObs = zeros(1,nAssims);
% for kk=1:nAssims
%     z(:,kk) = H*y(:,kk*Gap+1);
%     R(:,kk) = (.2*abs(z(:,kk))).^2;
%     for ll=1:3
%         if R(ll,kk) < MinNoise(ll)^2
%             R(ll,kk) = MinNoise(ll)^2;
%         end
%     end
%     z(:,kk) = z(:,kk)+sqrt(R(:,kk)).*randn(length(z(:,kk)),1);
%     tObs(kk) = t(kk*Gap+1);
% end
% %% ------------------------------------------------------------------------
% 
% 
% FileName = strcat('Stoch_SetUp_Gap_',num2str(Gap));
% save(FileName)


clear
close all
clc
%% ------------------------------------------------------------------------


%% Dan's set up
%% ------------------------------------------------------------------------
xo = [0.4424 0.3260 0.0308]';

dt = .05;
Gap = 12;
dT = Gap*dt;
nAssims = 1500;
T = nAssims*Gap*dt;
sqrtQ = 1*sqrt(dt)*eye(3);

%% Background covariance
xo = [-2.7295;-0.6538;24.2952];
[yC,~]=RunModel(xo,dt,10000,sqrtQ);
B = cov(yC');
sqrtB = sqrtm(B);

[y,t]=RunModel(xo,dt,T,sqrtQ);

H = [1 0 0; 0 0 1];
nObs = size(H,1);

R = zeros(nObs,nAssims);
z = zeros(nObs,nAssims);
tObs = zeros(1,nAssims);
for kk=1:nAssims
    z(:,kk) = H*y(:,kk*Gap+1);
    R(:,kk) = [0.1;0.1];
    z(:,kk) = z(:,kk)+sqrt(R(:,kk)).*randn(length(z(:,kk)),1);
    tObs(kk) = t(kk*Gap+1);
end
FileName = strcat('Dan_SetUp_Stoch_Gap_',num2str(Gap));
save(FileName)
disp('done')
%% ------------------------------------------------------------------------

