function [Xrs,inds] = resamplingDWD(weights,X,P,n)
c = zeros(P,1);
inds = zeros(P,1);
for jj=2:P+1
    c(jj)=c(jj-1)+weights(jj-1);
end, clear jj

% sample it and get the stronger particle more often
ii=1; % initialize
Xrs=zeros(n,P);
u1=rand/P; % intialize
for jj=1:P
    u = u1+(jj-1)/P;
    while u>c(ii)
        ii=ii+1;
    end
    Xrs(:,jj) = X(:,ii-1);
    inds(jj) = ii-1;
end
