function C = getCov(n,L)
% squared exponential
C = zeros(n);
for ii=1:n
    for jj=ii:n
        dist = min(abs(ii-jj));
        C(ii,jj) = exp(-(dist)^2/(2*L^2));
    end
end
C = (C+C')-diag(diag(C));


