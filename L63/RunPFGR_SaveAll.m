function RunPFGR_SaveAll(Ne,infl,SetUp,Gap)

load(strcat('./SetUps/SetUp_',num2str(SetUp),'_Gap_',num2str(Gap),'.mat'))
   
xo = y(:,1);
sqrtPo = diag([.1 .1 .05]);

[Xa,xAll,rho,traceP] = PF_GR_saveAll(Ne,infl,z,xo,sqrtPo,dt,dT,sqrtQ,H,R);

rmse = sqrt(sum((Xa-y(:,Gap+1:Gap:end)).^2)/3);
rmse = mean(rmse(200:end)); % time average
spread = sqrt(traceP/3);
spread = mean(spread(200:end)); % time average

Filename = strcat('./Results/PFGRSaveAll_Results_SetUp_',num2str(SetUp),'_Gap_',num2str(Gap), ...
    '_Ne_',num2str(Ne),'_infl_',num2str(infl),'.mat');
save(Filename,'rmse','spread','rho')

