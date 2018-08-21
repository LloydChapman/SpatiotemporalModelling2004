function [K,K0]=KnlSprtParas(d,alpha,typ,n,i1,i2,i3,j1,j2,j3)

if strcmp(typ,'Cauchy')
    K=1./(1+(d/alpha).^2);
elseif strcmp(typ,'Exp')
    K=exp(-d/alpha);
elseif strcmp(typ,'Const')
    K=ones(n);
end

K(1:n+1:end)=0; % set diagonal (same individual) entries to 0
% Zero different para entries
K(i1,j1)=0;
K(i2,j2)=0;
K(i3,j3)=0;
K0=n/sum(K(:));
K=K0*K; % normalise kernel
