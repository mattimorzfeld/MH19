%% define spatial domain
xS = -25*pi;
xF =  25*pi;
%% define time domain
dt = .1;
Nt = tF/dt + 1;

%% model parameters
md = -1;
mn = -1;
Np = 1;
Ng = -1;
Um = 1;
Uz = -1;
as = 0.0005;

% spatial grid 
Nx = 128;
dx = (xF-xS)/Nx;
x = (xS:dx:xF-dx)';

%% Obs
r = .1;
Nobs = Nx/ob_step;
H = zeros(Nobs,Nx);
for j = 0:Nobs-1
    H(j+1,ob_step*j+1) = 1;
end


load(['high_resolution_climatology_',num2str(Nx),'.mat']) 
ZH  = ZH(:,floor(length(ZH)/2):end);
