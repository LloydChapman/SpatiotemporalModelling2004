function [K,K0]=Knl(d,alpha,typ,n)

if strcmp(typ,'Cauchy')
    K=1./(1+(d/alpha).^2);
elseif strcmp(typ,'Exp')
    K=exp(-d/alpha);
elseif strcmp(typ,'Const')
    K=ones(n);
end

K(1:n+1:end)=0; % set diagonal (same individual) entries to 0
K0=n/sum(K(:));
K=K0*K; % normalise kernel

