function d_half=HalfRskDstnce(p,K0,typ)

beta=p(:,1);
alpha=p(:,2);
epsilon=p(:,3);
delta=p(:,4);

if strcmp(typ,'Const')
    d_half=NaN(size(p,1),1);
elseif strcmp(typ,'Cauchy')
    d_half=alpha.*sqrt(-beta.*K0./(log((1+exp(-(beta.*K0+delta+epsilon)))/2)+epsilon)-1);
elseif strcmp(typ,'Exp')
    d_half=alpha.*log(-beta.*K0./(log((1+exp(-(beta.*K0+delta+epsilon)))/2)+epsilon));
end
d_half(any(imag(d_half),2))=NaN;
d_half(d_half<0)=0;