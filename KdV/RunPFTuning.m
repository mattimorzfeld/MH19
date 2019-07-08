function RunPFTuning(Ne,tF,ob_step,infl)

FileName = strcat('./ProblemSetUps/SynthExp_Gap_',num2str(tF), ...
                    '_obStep_',num2str(ob_step),'.mat');
load(FileName)

locradAll = 1.5:.5:3;

for aa=1:length(locradAll)
    locrad = locradAll(aa);
    X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
    [xPF,tracePPF,Xa,rho]=PF(X,Z,locrad,infl,H,R,Nx,nAssims, ...
        Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
    rmse = sqrt(sum((xPF-Xt).^2)/Nx);
    spread = sqrt(tracePPF/Nx);
    
    rmse = mean(rmse(20:end));
    spread = mean(spread(20:end));
    
    FileName = strcat('./Results/PFResults_Gap_',num2str(tF), ...
        '_obStep_',num2str(ob_step),...
        '_Ne_',num2str(Ne), ...
        '_infl_',num2str(infl), ...
        '_loc_',num2str(locrad),...
        '.mat');
    save(FileName,'rmse','spread','rho')
end