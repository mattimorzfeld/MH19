%%
clear
close all
clc
colors

ob_step = 2; 
nAssims = 1200;
for tF = [2 4 6 8] 
    FileName = strcat('./ProblemSetUps/SpinUp_Gap_',num2str(tF), ...
                    '_obStep_',num2str(ob_step),'.mat'); 
    load(FileName)
    ProblemSetUp
    %% True State and obs
    [Z,Xt,R] = SynthExp(Xt(:,end),r,H,Nobs,nAssims,Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
    
    %% Background
    mub = xEnKF(:,end);
    dsm = 1; % Localization tuning
    build_localization
    B = 1.*CL.*cov(Xa');
    Lb = chol(B)';
    
    FileName = strcat('./LongProblemSetups/SynthExp_Gap_',num2str(tF), ...
                    '_obStep_',num2str(ob_step),'.mat'); 
    save(FileName)            
end
