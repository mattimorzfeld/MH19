%%
clear
% close all
clc

% GapAll = [2 4 5 8 10 11 12];
% Ne = 1000;
% rmse_OPF = zeros(length(GapAll),1);
% spread_OPF = zeros(length(GapAll),1);
% rho_OPF = zeros(length(GapAll),1);
% 
% for yy=1:length(GapAll)
%     Gap = GapAll(yy);
%     fprintf('Gap: %g\n',Gap)
% %     load(strcat('Stoch_SetUp_Gap_',num2str(Gap),'.mat'))
%     load(strcat('Dan_SetUp_Stoch_Gap_',num2str(Gap),'.mat'))
%     xo = y(:,1);
%     sqrtPo = diag([.1 .1 .05]);
%     [XaOPF,xAllOPF,rhoOPF,traceP_OPF] = OPF(Ne,z,xo,sqrtPo,dt,dT,sqrtQ,H,R);
%     
%     MSE = sqrt(sum((XaOPF-y(:,Gap+1:Gap:end)).^2)/3);
%     SPREAD = sqrt(traceP_OPF/3);
%     
%     rmse_OPF(yy) = mean(MSE(500:end));
%     spread_OPF(yy) = mean(SPREAD(500:end));
%     rho_OPF(yy) = mean(rhoOPF(500:end));
%     fprintf('RMSE, spread, rho: %g, %g, %g\n',rmse_OPF(yy),spread_OPF(yy),rho_OPF(yy))
% end
% 
% save TunedOPF_Dan_Stoch.mat

% load TunedOPF_Stoch.mat
load TunedOPF_Dan_Stoch.mat
colors

figure(1)
hold on, plot(GapAll*dt,rmse_OPF,'.-','Color',Color(:,6),'LineWidth',2,'MarkerSize',25)
hold on, plot(GapAll*dt,spread_OPF,'.--','Color',Color(:,6),'LineWidth',2,'MarkerSize',25)
xlabel('Time interval between obs.')
ylabel('RMSE and spread')
set(gca,'FontSize',16)
set(gcf,'Color','w')
box off

figure(3)
hold on, plot(GapAll*dt,rho_OPF,'.-','Color',Color(:,6),'LineWidth',2,'MarkerSize',25)
xlabel('Time interval between obs.')
ylabel('\rho')
set(gca,'FontSize',16)
set(gcf,'Color','w')
box off
