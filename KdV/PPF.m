function [xPF,tracePPF]=PPF(Ne,a,r,xpf,Z,H,R,Nx,nAssims,Nt,md,mn,Um,Uz,as,Np,Ng,xS,xF,dt)
xPF = zeros(Nx,nAssims);
tracePPF = zeros(1,nAssims); 
X = zeros(Nx,Ne);

for kk=1:nAssims
    fprintf('Jonathans PF assim %g / %g \r',kk,nAssims)
    % forecast
    for jj=1:Ne
        Cv = run_model_RK3(xpf{jj}, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
        tmp = ifft( Cv.','Symmetric' );
        xpf{jj} = tmp(:,end)';
    end
    [~,xpf,eflag] = pf_update_MWR16_beta(xpf,Nx,Ne,H,Z(:,kk)',r,a,1:size(H,1),R(:,kk),0);
    if eflag == 0
        fprintf('Error in localized PF.\n')
    end

    %% save results
    for oo=1:Ne
        X(:,oo) = xpf{oo}';
    end  
    xPF(:,kk) = mean(X,2);
    tracePPF(kk) = trace(cov(X'));
end