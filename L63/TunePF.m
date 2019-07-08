%%
clear
% close all
clc

% % GapAll = [5 8 10 11 12 13 14 15];
% GapAll = [2 4 5 8 10 11 12];
% Ne = 1000;
% rmse_PF = zeros(length(GapAll),1);
% spread_PF = zeros(length(GapAll),1);
% rho_PF = zeros(length(GapAll),1);
% 
% for yy=1:length(GapAll)
%     Gap = GapAll(yy);
%     fprintf('Gap: %g\n',Gap)
% %     load(strcat('Stoch_SetUp_Gap_',num2str(Gap),'.mat'))
%     load(strcat('Dan_SetUp_Stoch_Gap_',num2str(Gap),'.mat'))
%     xo = y(:,1);
%     sqrtPo = diag([.1 .1 .05]);
%     [XaPF,xAllPF,rhoPF,traceP_PF] = PF(Ne,z,xo,sqrtPo,dt,dT,sqrtQ,H,R);
%     
%     MSE = sqrt(sum((XaPF-y(:,Gap+1:Gap:end)).^2)/3);
%     SPREAD = sqrt(traceP_PF/3);
%     
%     rmse_PF(yy) = mean(MSE(500:end));
%     spread_PF(yy) = mean(SPREAD(500:end));
%     rho_PF(yy) = mean(rhoPF(500:end));
%     clc
% end
% 
% save TunedPF_Dan_Stoch.mat

% load TunedPF_Stoch.mat
load TunedPF_Dan_Stoch.mat
colors

figure(1)
hold on, plot(GapAll*dt,rmse_PF,'.-','Color',Color(:,5),'LineWidth',2,'MarkerSize',25)
hold on, plot(GapAll*dt,spread_PF,'.--','Color',Color(:,5),'LineWidth',2,'MarkerSize',25)
xlabel('Time interval between obs.')
ylabel('RMSE and spread')
set(gca,'FontSize',16)
set(gcf,'Color','w')
box off

figure(3)
hold on, plot(GapAll*dt,rho_PF,'.-','Color',Color(:,5),'LineWidth',2,'MarkerSize',25)
xlabel('Time interval between obs.')
ylabel('\rho')
set(gca,'FontSize',16)
set(gcf,'Color','w')
box off
