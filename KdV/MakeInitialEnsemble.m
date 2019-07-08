function iEns = MakeInitialEnsemble(mu,Lb,Ne,Nx)

iEns = zeros(Nx,Ne);
for kk=1:Ne
    go = 1;
    while go == 1
        xi = Lb*randn(Nx,1);
        sp = mu + xi;
        sm = mu - xi;
        if min(sp)>-1
            iEns(:,kk) = sp;
            go = 0;
        elseif min(sm)>-1
            iEns(:,kk) = sm;
            go = 0;
        end
    end
end