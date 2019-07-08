%% create sinusoid matrix (basis functions)
% ... because we use sinusoids our correlation matrix will be periodic!

E = zeros(Nx);
zj = (2*pi/Nx:2*pi/Nx:2*pi)';
e0 = ones(Nx,1)./sqrt(Nx);
E(:,1) = e0;
for i = 1:Nx/2-1,
  tt = cos(i*zj);
  tn = sqrt(tt'*tt);
  E(:,2*i) = tt./tn;
  tt = sin(i*zj);
  tn = sqrt(tt'*tt);
  E(:,2*i+1) = tt./tn;
end
tt = cos(Nx*zj/2);
tn = sqrt(tt'*tt);
E(:,Nx) = tt./tn;

%% create "eigenvalues"

g = ones(1,Nx);
b = 1./dsm^2;
for i = 1:(Nx-2)/2
  g(2*i)   = exp(-b*i^2);
  g(2*i+1) = g(2*i);
end
g(Nx) = exp(-b*(Nx/2)^2);

a = Nx/sum(g);
g = a*g;
G = diag(g);

%% create localization matrix

CL = E*G*E';

% % ensure distant correlations are identically zero
% for j = 1:Nx
%     for p = 1:Nx
%         if abs(CL(j,p)) < 1e-4
%             CL(j,p) = 0;
%         end
%     end
% end