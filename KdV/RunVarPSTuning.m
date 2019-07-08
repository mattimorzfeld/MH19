function RunVarPSTuning(Ne,tF,ob_step,infl,OuterLoops)

FileName = strcat('./ProblemSetUps/SynthExp_Gap_',num2str(tF), ...
                    '_obStep_',num2str(ob_step),'.mat');
load(FileName)

locradAll = 1.:.5:3;

for aa=1:length(locradAll)
    locrad = locradAll(aa);
    [xvarPSww,tracePvarPSww,~]=varPSww(Ne,mub,Lb,Z,locrad,infl,H,R,Nx,nAssims, ...
        Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt, OuterLoops);
    rmse = sqrt(sum((xvarPSww-Xt).^2)/Nx);
    spread = sqrt(tracePvarPSww/Nx);
    
    rmse = mean(rmse(20:end));
    spread = mean(spread(20:end));
    
    FileName = strcat('./Results/varPSResults_Gap_',num2str(tF), ...
        '_obStep_',num2str(ob_step),...
        '_Ne_',num2str(Ne), ...
        '_infl_',num2str(infl), ...
        '_loc_',num2str(locrad),...
        '_OuterLoops',num2str(OuterLoops),...
        '.mat');
    save(FileName,'rmse','spread')
end
  
