%%
clear
close all
clc
colors

ob_step = 2; 
nAssims = 120;

for tF = [2 4 6 8]
    ProblemSetUp
    %% True State and obs
    [Z,Xt,R] = SynthExp(ZH(:,end),r,H,Nobs,nAssims,Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
    
    %% Background
    mu = Xt(:,1);
    B = cov(ZH');
    Lb = chol(B)';
    
    %% EnKF
    Ne = 50;
    locrad = 3;
    infl = 1.2;
    X = MakeInitialEnsemble(mu,Lb,Ne,Nx);
    
    [xEnKF,tracePEnKF,Xa]=EnKF(X,Z,locrad,infl,H,R,Nx,Nobs,nAssims, ...
        Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
    rmse_EnKF = sqrt(sum((xEnKF-Xt).^2)/Nx);
    spread_EnKF = sqrt(tracePEnKF/Nx);
    
    %%
    figure
    plot(rmse_EnKF)
    hold on, plot(spread_EnKF,'--')
    
    %% animation
%     figure
%     set(gcf,'Color','w')
%     for kk=1:nAssims
%         plot(x,Xt(:,kk),'Color',Color(:,1),'LineWidth',2);
%         hold on, plot(x,xEnKF(:,kk),'Color',Color(:,5),'LineWidth',2);
%         hold on,plot(x(1:ob_step:end),Z(:,kk),'.','Color',Color(:,4),'MarkerSize',20);
%         axis([-80 80 -2 6])
%         box off
%         set(gca,'FontSize',12)
%         drawnow
%         hold off
%     end
    
    FileName = strcat('./ProblemSetUps/SpinUp_Gap_',num2str(tF), ...
        '_obStep_',num2str(ob_step),'.mat');
    save(FileName,'Xa', 'Xt', 'xEnKF','ob_step','tF')
end
