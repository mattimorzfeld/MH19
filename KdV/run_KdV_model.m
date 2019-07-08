clearvars
close all
clc

format long

%% define spatial domain

xS = -25*pi;
xF =  25*pi;

%% length of forecast

tF = 200;
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
Nx = 256;
dx = (xF-xS)/Nx;
x = (xS:dx:xF-dx)';

%% make initial condition

IC = zeros(Nx,1);
for j = 1:Nx
    IC(j) = exp(-1/2/100*x(j)^2 );
end

%% run model
Cv = run_model_RK3(IC, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt);
Xt = ifft( Cv.','Symmetric' );

figure(1)
clf
contourf(x,t,Xt',21,'LineColor','none')
colorbar