clearvars
format long

%% define spatial domain

xS = -25*pi;
xF =  25*pi;

%% length of forecast

tF = 1e3;
dt = 1e-2;
Nt = tF/dt + 1;
t = 0:dt:tF;

%% model parameters

md = -1;
mn = -1;
Np = 1;
Ng = -1;
Um = 1;
Uz = -1;
as = 0.0005;

%% spatial grid
Nx = 128;
dx = (xF-xS)/Nx;
x = (xS:dx:xF-dx)';

%% make initial condition

IC = zeros(Nx,1);
for j = 1:Nx
    IC(j) = exp(-1/2/100*x(j)^2 );
end

%% run model
Cv = run_model_RK3(IC, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
ZH = ifft( Cv.','Symmetric' );

save(['high_resolution_climatology_',num2str(Nx),'.mat'],'ZH') 


%%
ind = 80000;
ZHn = ZH(:,ind:end);
tn= t(:,ind:end);
[X,T]=meshgrid(x,tn);
mesh(X,T,ZHn')