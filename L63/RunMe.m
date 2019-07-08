%%
clear
close all
clc
%% ------------------------------------------------------------------------

SetUp = 2;
Gap = 2;

load(strcat('./SetUps/SetUp_',num2str(SetUp),'_Gap_',num2str(Gap),'.mat'))


%% EnKF (stochastic)
%% ------------------------------------------------------------------------
if 1
Ne = 50;
infl = 1.2;
[XaEnKF,xAllEnKF,traceP_EnKF] = EnKF(Ne,infl,z,dt,dT,sqrtQ,H,R);

rmse_EnKF = sqrt(sum((XaEnKF-y(:,Gap+1:Gap:end)).^2)/3);
spread_EnKF = sqrt(traceP_EnKF/3);

MakePlot(y,z,t,tObs,xAllEnKF,rmse_EnKF,spread_EnKF,0,'EnKF')
end
%% ------------------------------------------------------------------------


%% varPS 
%% ------------------------------------------------------------------------
if 0
Ne = 100;
mu = xo;
infl = 1.9;
[Xa_varPS,xAll_varPS,rho_varPS,traceP_varPS] = varPS(Ne,z,mu,sqrtB,infl,dt,dT,H,R);

rmse_varPS = sqrt(sum((Xa_varPS-y(:,Gap+1:Gap:end)).^2)/3);
spread_varPS = sqrt(traceP_varPS/3);

MakePlot(y,z,t,tObs,xAll_varPS,rmse_varPS,spread_varPS,rho_varPS,'varPS')
end
%% ------------------------------------------------------------------------


%% varPS without weights
%% ------------------------------------------------------------------------
if 0
Ne = 100;
mu = xo;
infl = 1.2;
[Xa_varPSnw,xAll_varPSnw,traceP_varPSnw] = varPSnw(Ne,z,mu,sqrtB,infl,dt,dT,H,R);

rmse_varPSnw = sqrt(sum((Xa_varPSnw-y(:,Gap+1:Gap:end)).^2)/3);
spread_varPSnw = sqrt(traceP_varPSnw/3);

MakePlot(y,z,t,tObs,xAll_varPSnw,rmse_varPSnw,spread_varPSnw,0,'varPS no weights')
end
%% ------------------------------------------------------------------------


%% EDA
%% ------------------------------------------------------------------------
if 0
Ne = 20;
mu = xo;
infl = 1.2;
[Xa_EDA,xAll_EDA,traceP_EDA] = EDA(Ne,z,mu,sqrtB,infl,dt,dT,H,R);

rmse_EDA = sqrt(sum((Xa_EDA-y(:,Gap+1:Gap:end)).^2)/3);
spread_EDA = sqrt(traceP_EDA/3);

MakePlot(y,z,t,tObs,xAll_EDA,rmse_EDA,spread_EDA,0,'EDA')
end
%% ------------------------------------------------------------------------

%% Standard PF
%% ------------------------------------------------------------------------
if 0
Ne = 1000;
xo = y(:,1);
sqrtPo = diag([.1 .1 .05]);
[XaPF,xAllPF,rhoPF,traceP_PF] = PF(Ne,z,xo,sqrtPo,dt,dT,sqrtQ,H,R);

rmse_PF = sqrt(sum((XaPF-y(:,Gap+1:Gap:end)).^2)/3);
spread_PF =sqrt(traceP_PF/3);
%%
MakePlot(y,z,t,tObs,xAllPF,rmse_PF,spread_PF,rhoPF,'PF')
end
%% ------------------------------------------------------------------------


%% Standard PF with OPF step
%% ONLY WORKS WITH Q =\= 0
%% ------------------------------------------------------------------------
if 0
Ne = 200;
xo = y(:,1);
sqrtPo = diag([.1 .1 .05]);
[XaOPF,xAllOPF,rhoOPF,traceP_OPF] = OPF(Ne,z,xo,sqrtPo,dt,dT,sqrtQ,H,R);

rmse_OPF = sqrt(sum((XaOPF-y(:,Gap+1:Gap:end)).^2)/3);
spread_OPF = sqrt(traceP_OPF/3);

MakePlot(y,z,t,tObs,xAllOPF,rmse_OPF,spread_OPF,rhoOPF,'OPF')
end
%% ------------------------------------------------------------------------


%% Standard PF with Gaussian replenishing
%% ------------------------------------------------------------------------
if 0
Ne = 500;
xo = y(:,1);
sqrtPo = diag([.1 .1 .05]);
infl = 1.2;
[XaPF_GR,xAllPF_GR,rhoPF_GR,traceP_PF_GR] = PF_GR(Ne,infl,z,xo,sqrtPo,dt,dT,sqrtQ,H,R);

rmse_PF_GR = sqrt(sum((XaPF_GR-y(:,Gap+1:Gap:end)).^2)/3);
spread_PF_GR = sqrt(traceP_PF_GR/3);

MakePlot(y,z,t,tObs,xAllPF_GR,rmse_PF_GR,spread_PF_GR,rhoPF_GR,'PF w GR')
end
%% ------------------------------------------------------------------------

%% DWD PF
%% ------------------------------------------------------------------------
if 1
Ne = 500;
xo = y(:,1);
sqrtPo = diag([.1 .1 .05]);
[XaPF_DWD,xAllPF_DWD,rhoPF_DWD,traceP_PF_DWD] = PF_DWD(Ne,1e-1,z,xo,sqrtPo,dt,dT,sqrtQ,H,R);

rmse_PF_DWD= sqrt(sum((XaPF_DWD-y(:,Gap+1:Gap:end)).^2)/3);
spread_PF_DWD = sqrt(traceP_PF_DWD/3);

MakePlot(y,z,t,tObs,xAllPF_DWD,rmse_PF_DWD,spread_PF_DWD,rhoPF_DWD,'DWD PF')
end
%% ------------------------------------------------------------------------



%% Display results
%% ------------------------------------------------------------------------
disp(' '),disp(' '),disp(' ')

% fprintf('PF - rmse / spead: %g / %g \n',mean(rmse_PF),mean(sqrt(traceP_PF/3)))
% fprintf('OPF - rmse / spead: %g / %g \n',mean(rmse_OPF),mean(sqrt(traceP_OPF/3)))
% fprintf('PF GR - rmse / spead: %g / %g \n',mean(rmse_PF_GR),mean(sqrt(traceP_PF_GR/3)))
fprintf('DWD PF - rmse / spead: %g / %g \n',mean(rmse_PF_DWD),mean(sqrt(traceP_PF_DWD/3)))

% fprintf('varPS - rmse / spead: %g / %g \n',mean(rmse_varPS),mean(sqrt(traceP_varPS/3)))
% fprintf('varPS no weights - rmse / spead: %g / %g \n',mean(rmse_varPSnw),mean(sqrt(traceP_varPSnw/3)))

% fprintf('EDA - rmse / spead: %g / %g \n',mean(rmse_EDA),mean(sqrt(traceP_EDA/3)))

fprintf('EnKF - rmse / spead: %g / %g \n',mean(rmse_EnKF),mean(sqrt(traceP_EnKF/3)))
% ------------------------------------------------------------------------

