%%
clear 
clc
colors


SetUp = 1;
Gap = [14];

% close all
% Method = 'EnKF';
% ColorNumber = 2;
% infl = .9:.05:2;
% NeAll = [20 50];% 100 200 500];

Method = 'PFDWD';
ColorNumber = 4;
infl = .9:.05:2;
NeAll = [100 200];% 100 200 500];



for kk=1:length(NeAll)
    Ne = NeAll(kk);
    RMSEAll = zeros(length(infl),1);
    SpreadAll = zeros(length(infl),1);
    for ll = 1:length(infl)
        Filename = strcat('./Results/',Method,'_Results_SetUp_',num2str(SetUp),'_Gap_',num2str(Gap), ...
            '_Ne_',num2str(Ne),'_infl_',num2str(infl(ll)),'.mat');
        load(Filename)
        RMSEAll(ll) = mean(rmse(200:end));
        SpreadAll(ll) = mean(spread(200:end));
    end
    [~,ind] = min(RMSEAll);
    tuned_RMSE(kk) = RMSEAll(ind);
    tuned_Spread(kk) = SpreadAll(ind);
end
figure(1)
plot(NeAll,tuned_RMSE,'.-','Color',Color(:,ColorNumber),'LineWidth',2,'MarkerSize',30)
hold on,plot(NeAll,tuned_Spread,'.','Color',Color(:,ColorNumber),'MarkerSize',30)
hold on,plot(NeAll,tuned_Spread,'--','Color',Color(:,ColorNumber),'LineWidth',2)
xlabel('Time interval between obs.')
ylabel('RMSE and spread')
set(gcf,'Color','w')
set(gca,'FontSize',16)
box off


