%%
clear
close all
clc
colors

tF = 8;
ob_step = 2;
FileName = strcat('./ProblemSetUps/SynthExp_Gap_',num2str(tF), ...
                    '_obStep_',num2str(ob_step),'.mat');
load(FileName)


%% Jonathan'sPF
%% ----------------------------------------------------------------------
% Ne = 50;
% r = 2;
% a = 0.1;
% X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
% xpf = cell(Ne,1);
% for kk=1:Ne
%     xpf{kk} = X(:,kk);
% end
%     
% [xJPPF,tracePJPPF]=PPF(Ne,a*Ne,r,xpf,Z,H,R,Nx,nAssims,Nt,md,mn,Um,Uz,as,Np,Ng,xS,xF,dt);
% rmse_JPPF = sqrt(sum((xJPPF-Xt).^2)/Nx);
% spread_JPPF = sqrt(tracePJPPF/Nx);
% %%
% figure(1)
% hold on, plot(rmse_JPPF,'Color',Color(:,3),'LineWidth',2)
% hold on, plot(spread_JPPF,'--','Color',Color(:,3),'LineWidth',2)
% drawnow
%% ----------------------------------------------------------------------


%% PF
%% ----------------------------------------------------------------------
% Ne = 1000;
% locrad = 2;
% infl = 1.2;
% X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
% [xPF,tracePPF,Xa,rho]=PF(X,Z,locrad,infl,H,R,Nx,nAssims, ...
%                                     Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
% rmse_PF = sqrt(sum((xPF-Xt).^2)/Nx);
% spread_PF = sqrt(tracePPF/Nx);
% %%
% figure(1)
% hold on, plot(rmse_PF,'Color',Color(:,4),'LineWidth',2)
% hold on, plot(spread_PF,'--','Color',Color(:,4),'LineWidth',2)
% drawnow
%% ----------------------------------------------------------------------

%% varPS without weights
%% ----------------------------------------------------------------------
% Ne = 20;
% locrad = 2;
% infl = 1.2;
% OuterLoops = 2;
% [xvarPSww,tracePvarPSww,Xa_varPSww]=varPSww(Ne,mub,Lb,Z,locrad,infl,H,R,Nx,nAssims, ...
%                                     Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt, OuterLoops);
% rmse_varPSww = sqrt(sum((xvarPSww-Xt).^2)/Nx);
% spread_varPSww = sqrt(tracePvarPSww/Nx);
% %%
% figure(1)
% hold on, plot(rmse_varPSww,'Color',Color(:,1),'LineWidth',2)
% hold on, plot(spread_varPSww,'--','Color',Color(:,1),'LineWidth',2)
% drawnow
%% ----------------------------------------------------------------------

%% localized PF
%% ----------------------------------------------------------------------
% Ne = 100;
% locrad = 2;
% infl = 1.1;
% HowManyObs = 3;
% X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
% [xlPF,tracePlPF,Xa]=lPF(X,Z,locrad,infl,H,R,Nx,nAssims, ...
%                                     Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt,ob_step,HowManyObs);
% rmse_lPF = sqrt(sum((xlPF-Xt).^2)/Nx);
% spread_lPF = sqrt(tracePlPF/Nx);
% %%
% figure(1)
% hold on, plot(rmse_lPF,'Color',Color(:,4),'LineWidth',2)
% hold on, plot(spread_lPF,'--','Color',Color(:,4),'LineWidth',2)
% drawnow
%% ----------------------------------------------------------------------


%% EDA
%% ----------------------------------------------------------------------
% Ne = 20;
% locrad = 2;
% infl = 1.2;
% OuterLoops = 0;
% [xEDA,tracePEDA]=EDA(Ne,mub,Lb,Z,locrad,infl,H,R,Nx,nAssims, ...
%                                     Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt, OuterLoops);
% rmse_EDA = sqrt(sum((xEDA-Xt).^2)/Nx);
% spread_EDA = sqrt(tracePEDA/Nx);
% %
% figure(1)
% hold on, plot(rmse_EDA,'Color',Color(:,6),'LineWidth',2)
% hold on, plot(spread_EDA,'--','Color',Color(:,6),'LineWidth',2)
% drawnow
%% ----------------------------------------------------------------------


%% EnKF
%% ----------------------------------------------------------------------
Ne = 20;
locrad = 2;
infl = 1.2;
X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
[xEnKF,tracePEnKF,Xa]=EnKF(X,Z,locrad,infl,H,R,Nx,Nobs,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
rmse_EnKF = sqrt(sum((xEnKF-Xt).^2)/Nx);
spread_EnKF = sqrt(tracePEnKF/Nx);

%
figure(1)
hold on, plot(rmse_EnKF,'Color',Color(:,2),'LineWidth',2)
hold on, plot(spread_EnKF,'--','Color',Color(:,2),'LineWidth',2)
drawnow
%% ----------------------------------------------------------------------

%% EnKS
%% ----------------------------------------------------------------------
Ne = 20;
locrad = 2;
infl = 1.2;
X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
[xEnKS,tracePEnKS,Xa]=EnKS(X,Z,locrad,infl,H,R,Nx,Nobs,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
rmse_EnKS = sqrt(sum((xEnKS-Xt).^2)/Nx);
spread_EnKS = sqrt(tracePEnKS/Nx);

%%
figure(1)
hold on, plot(rmse_EnKS,'Color',Color(:,6),'LineWidth',2)
hold on, plot(spread_EnKS,'--','Color',Color(:,6),'LineWidth',2)
drawnow
%% ----------------------------------------------------------------------


%% Print results
%% ----------------------------------------------------------------------
% fprintf('localized PF RMSE = %g,\n',mean(rmse_lPF(20:end)))
fprintf('EnKF RMSE = %g,\n',mean(rmse_EnKF(20:end)))
fprintf('EnKS RMSE = %g,\n',mean(rmse_EnKS(20:end)))
% fprintf('Loc.PF RMSE = %g,\n',mean(rmse_JPPF(20:end)))
% fprintf('varPS RMSE = %g,\n',mean(rmse_varPSww(20:end)))
% fprintf('EDA RMSE = %g,\n',mean(rmse_EDA(20:end)))

%% ----------------------------------------------------------------------




%% animation
%% ----------------------------------------------------------------------
% figure
% set(gcf,'Color','w')
% for kk=1:nAssims
%     plot(x,Xt(:,kk),'Color',Color(:,1),'LineWidth',2);
%     hold on, plot(x,xEnKF(:,kk),'Color',Color(:,5),'LineWidth',2);
%     hold on,plot(x(1:ob_step:end),Z(:,kk),'.','Color',Color(:,4),'MarkerSize',20);
%     axis([-80 80 -2 6])
%     box off
%     set(gca,'FontSize',12)
%     drawnow
%     hold off
% end

