%%
clear 
clc
colors

% close all

Method = 'EnKF';
ColorNumber = 2;
infl = .9:.05:2;
Ne = 20;

% Method = 'PFDWD';
% ColorNumber = 1;
% infl = .9:.05:2;
% Ne = 200;

% Method = 'PFGR';
% ColorNumber = 3;
% infl = .9:.05:2;
% Ne = 100;

% Method = 'varPSnw';
% ColorNumber = 4;
% infl = .9:.05:2;
% Ne = 20;

% Method = 'varPS';
% ColorNumber = 5;
% infl = .9:.05:2;
% Ne = 20;

SetUps = 1:4;
Gaps = [2 4 6 8 10 12 14];

for jj = 1:length(SetUps)
    SetUp = SetUps(jj);
    
    tuned_RMSE_AtGap = zeros(length(Gaps),1);
    tuned_Spread_AtGap = zeros(length(Gaps),1);
    dTs = zeros(length(Gaps),1);
    for kk = 1:length(Gaps)
        Gap = Gaps(kk);
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
        tuned_RMSE_AtGap(kk) = RMSEAll(ind);
        tuned_Spread_AtGap(kk) = SpreadAll(ind);
        dTs(kk) = Gap*.05;
    end
    
    figure(jj)
    plot(dTs,tuned_RMSE_AtGap,'.-','Color',Color(:,ColorNumber),'LineWidth',2,'MarkerSize',30)
    hold on,plot(dTs,tuned_Spread_AtGap,'.','Color',Color(:,ColorNumber),'MarkerSize',30)
    hold on,plot(dTs,tuned_Spread_AtGap,'--','Color',Color(:,ColorNumber),'LineWidth',2)
    xlabel('Time interval between obs.')
    ylabel('RMSE and spread')
    set(gcf,'Color','w')
    set(gca,'FontSize',16)
    box off
end