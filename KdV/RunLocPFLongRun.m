function RunLocPFLongRun(Ne,tF,ob_step,infl,locrad)

FileName = strcat('./LongProblemSetups/SynthExp_Gap_',num2str(tF), ...
                    '_obStep_',num2str(ob_step),'.mat');
load(FileName)
     
X = MakeInitialEnsemble(mub,Lb,Ne,Nx);
xpf = cell(Ne,1);
for kk=1:Ne
    xpf{kk} = X(:,kk);
end
[xJPPF,tracePJPPF]=PPF(Ne,infl*Ne,r,xpf,Z,H,R,Nx,nAssims,Nt,md,mn,Um,Uz,as,Np,Ng,xS,xF,dt);
rmse = sqrt(sum((xJPPF-Xt).^2)/Nx);
spread = sqrt(tracePJPPF/Nx);

rmse = mean(rmse(200:end));
spread = mean(spread(200:end));

FileName = strcat('./ResultsTuned/LocPFResults_Gap_',num2str(tF), ...
    '_obStep_',num2str(ob_step),...
    '_Ne_',num2str(Ne), ...
    '_infl_',num2str(infl), ...
    '_loc_',num2str(locrad),...
    '.mat');
save(FileName,'rmse','spread')
  