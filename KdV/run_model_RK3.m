function Cv = run_model_RK3(IC, Nx, Nt, md, mn, Um, Uz, as, Np, Ng, xS, xF, dt)

%% spatial domain

% grid
NxB = 2*Nx;
Ncor = NxB/Nx;

dx = (xF - xS)/(Nx - 0);
x  = (xS:dx:xF-dx)';
dxB = (xF-xS)/(NxB - 0);
xB = (xS:dxB:xF-dxB)';

mp_x = ( Um + Uz*exp( -as*xB.^2 ) )*Np;
mg_x = -2*Uz*as*Ng*xB.*exp( -as*xB.^2 );

%% sponge

sp = -2;
Lsp = 2*pi;

spongeB = sp*exp( -(1/Lsp^2)*( xB - xS ).^2 ) ...
        + sp*exp( -(1/Lsp^2)*( xB - xF ).^2 );

%% create Fourier factors for derivatives

for j = 1:Nx/2+1
  ddx(j) = 1i*2*pi*(j-1)/(xF-xS);
end

for j = Nx/2+2:Nx
  ddx(j) = -ddx(Nx-j+2);
end
        
%% intialize 

Cv = zeros(Nt,Nx);
Cnonlin_adv = zeros(1,Nx);

%% create initial state

Cv(1,:) = fft( IC ).';

%% integrate with RK3

for n = 1:Nt-1

    %% step 1
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = Cv(n,1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = Cv(n,Nx/2+2:Nx);
    AB = ifft(Cwrk,'symmetric')*Ncor;
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = ddx(1:Nx/2+1).*Cv(n,1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = ddx(Nx/2+2:Nx).*Cv(n,Nx/2+2:Nx);
    dAdxB = ifft(Cwrk,'symmetric')*Ncor;
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = ddx(1:Nx/2+1).^3.*Cv(n,1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = ddx(Nx/2+2:Nx).^3.*Cv(n,Nx/2+2:Nx);
    d3Adx3B = ifft(Cwrk,'symmetric')*Ncor;
        
    % sponge
    BC = spongeB.*AB;
    
    nonlin_adv = -( md*d3Adx3B + ( mp_x + mn*AB ).*dAdxB + mg_x.*AB ) + BC;
    Cwrk = zeros(NxB,1);
    Cwrk = fft(nonlin_adv);
    Cnonlin_adv(1:Nx/2+1) = Cwrk(1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nx/2+2:Nx) = Cwrk(NxB-Nx/2+2:NxB)/Ncor;
    
    % tendency

    C1_tendency = Cnonlin_adv;

    % 

    C1 = Cv(n,:) + dt*C1_tendency/3;
    
    %% step 2
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = C1(1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = C1(Nx/2+2:Nx);
    AB = ifft(Cwrk,'symmetric')*Ncor;
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = ddx(1:Nx/2+1).*C1(1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = ddx(Nx/2+2:Nx).*C1(Nx/2+2:Nx);
    dAdxB = ifft(Cwrk,'symmetric')*Ncor;
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = ddx(1:Nx/2+1).^3.*C1(1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = ddx(Nx/2+2:Nx).^3.*C1(Nx/2+2:Nx);
    d3Adx3B = ifft(Cwrk,'symmetric')*Ncor;
    
    % sponge
    BC = spongeB.*AB;
    
    nonlin_adv = -( md*d3Adx3B + ( mp_x + mn*AB ).*dAdxB + mg_x.*AB ) + BC;
    Cwrk = zeros(NxB,1);
    Cwrk = fft(nonlin_adv);
    Cnonlin_adv(1:Nx/2+1) = Cwrk(1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nx/2+2:Nx) = Cwrk(NxB-Nx/2+2:NxB)/Ncor;
    
    % tendency

    C2_tendency = Cnonlin_adv - 5*C1_tendency/9;

    %

    C2 = C1 + 15*dt*C2_tendency/16;
    
    %% step 3
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = C2(1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = C2(Nx/2+2:Nx);
    AB = ifft(Cwrk,'symmetric')*Ncor;
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = ddx(1:Nx/2+1).*C2(1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = ddx(Nx/2+2:Nx).*C2(Nx/2+2:Nx);
    dAdxB = ifft(Cwrk,'symmetric')*Ncor;
    Cwrk = zeros(NxB,1);
    Cwrk(1:Nx/2+1) = ddx(1:Nx/2+1).^3.*C2(1:Nx/2+1);
    Cwrk(NxB-Nx/2+2:NxB) = ddx(Nx/2+2:Nx).^3.*C2(Nx/2+2:Nx);
    d3Adx3B = ifft(Cwrk,'symmetric')*Ncor;
    
    % sponge
    BC = spongeB.*AB;
    
    nonlin_adv = -( md*d3Adx3B + ( mp_x + mn*AB ).*dAdxB + mg_x.*AB ) + BC;
    Cwrk = zeros(NxB,1);
    Cwrk = fft(nonlin_adv);
    Cnonlin_adv(1:Nx/2+1) = Cwrk(1:Nx/2+1)/Ncor;
    Cnonlin_adv(Nx/2+2:Nx) = Cwrk(NxB-Nx/2+2:NxB)/Ncor;
    
    % tendency

    C3_tendency = Cnonlin_adv - 153*C2_tendency/128;

    % 
    
    Cv(n+1,:) = C2 + 8*dt*C3_tendency/15;
    
end