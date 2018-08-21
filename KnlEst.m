function Ke=KnlEst(d,alpha,typ)

if strcmp(typ,'Cauchy')
    Ke=1./(1+(d./alpha).^2);
elseif strcmp(typ,'Exp')
    Ke=exp(-d./alpha);
elseif strcmp(typ,'Const')
    Ke=ones(1,numel(d));
end