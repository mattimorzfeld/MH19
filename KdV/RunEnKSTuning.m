function RunEnKSTuning(Ne,tF,ob_step)

FileName = strcat('./ProblemSetUps/SynthExp_Gap_',num2str(tF), ...
                    '_obStep_',num2str(ob_step),'.mat');
load(FileName)
inflAll = .9:.1:2;
locradAll = 1:.5:4;
for qq=1:length(inflAll)
    infl = inflAll(qq);
    for pp=1:length(locradAll)
        locrad = locradAll(pp);
        
        X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
        [xEnKS,tracePEnKS,~]=EnKS(X,Z,locrad,infl,H,R,Nx,Nobs,nAssims, ...
                                    Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        rmse = sqrt(sum((xEnKS-Xt).^2)/Nx);
        spread = sqrt(tracePEnKS/Nx);
        
        rmse = mean(rmse(20:end));
        spread = mean(spread(20:end));
        
        FileName = strcat('./Results/EnKSResults_Gap_',num2str(tF), ...
                                '_obStep_',num2str(ob_step),...
                                '_Ne_',num2str(Ne), ...
                                '_infl_',num2str(infl), ...
                                '_loc_',num2str(locrad),'.mat');
        save(FileName,'rmse','spread')
    end
end    