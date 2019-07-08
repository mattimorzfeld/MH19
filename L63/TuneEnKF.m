function TuneEnKF(Ne,inflAll,SetUp,Gap)

load(strcat('./SetUps/SetUp_',num2str(SetUp),'_Gap_',num2str(Gap),'.mat'))

for kk = 1:length(inflAll)
    infl = inflAll(kk);
    [Xa,xAll,traceP] = EnKF(Ne,infl,z,dt,dT,sqrtQ,H,R);
    
    rmse = sqrt(sum((Xa-y(:,Gap+1:Gap:end)).^2)/3);
    rmse = mean(rmse(200:end)); % time average
    spread = sqrt(traceP/3);
    spread = mean(spread(200:end)); % time average
    
    Filename = strcat('./Results/EnKF_Results_SetUp_',num2str(SetUp),'_Gap_',num2str(Gap), ...
        '_Ne_',num2str(Ne),'_infl_',num2str(infl),'.mat');
    save(Filename,'rmse','spread')
%     save(Filename,'y','z','t','tObs','xAll','rmse','spread','Xa')
end