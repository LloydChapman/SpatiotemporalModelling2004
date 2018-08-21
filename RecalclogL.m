function LL=RecalclogL(n,tmax,I,RL,tE,tI,tRorD,tRL,tRLR,preB,tDm,prevK,IP,r1,p1,h0,d,d0,beta,alpha,epsilon,delta,typ,tLm)
tEm=false(n,tmax);
tEm((tE-1)*n+I)=1;
if isempty(tLm)
    S=1-preB-max(cumsum(tEm,2),cumsum(tDm,2)); % don't remove LST+ individuals
else
    S=1-preB-max(max(cumsum(tEm,2),cumsum(tDm,2)),cumsum(tLm,2)); % remove LST+ individuals from susceptibles
end
S(prevK,:)=0; % remove previous KA cases from susceptibles

h=zeros(n,tmax);
for i=1:numel(I)
    h(I(i),tE(i):tI(I(i))-1)=h0;
    h(I(i),tI(I(i)):tRorD(I(i))-1)=1;
end

for i=1:numel(RL)
    h(RL(i),tRL(RL(i)):min(tRLR(RL(i))-1,tmax,'omitnan'))=1;
end

K=Knl(d,alpha,typ,n);
rate=beta*K+delta*d0;
lambda=rate(:,I)*h(I,:)+epsilon;

LL=logL(S,lambda,tEm,IP,r1,p1);