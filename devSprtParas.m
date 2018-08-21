function [D,LL]=devSprtParas(p,p1,z,n,tmax,NONR,RNO,ONR,tIsNONR,tRsNONR,tIsRNO,tRsONR,I,RL,tEs,tI,tRorD,tRL,tRLR,preB,tDm,prevK,IPm,r1,h0,d,d0,typ,tLm,ss,i1,i2,i3,j1,j2,j3)

if strcmp(ss,'mean')
    beta=mean(p(z,1));
    alpha=mean(p(z,2));
    epsilon=mean(p(z,3));
    delta=mean(p(z,4));
    p1=mean(p1(z));
elseif strcmp(ss,'mode')
    beta=MAP(p(z,1));
    alpha=MAP(p(z,2));
    epsilon=MAP(p(z,3));
    delta=MAP(p(z,4));
    p1=MAP(p1(z));
end

LL=zeros(numel(z),1);
for i=1:numel(z)
%     i
    k=z(i);
    tI(NONR)=tIsNONR(:,k);
    tRorD(NONR)=tRsNONR(:,k);
    tI(RNO)=tIsRNO(:,k);
    tRorD(ONR)=tRsONR(:,k);
    LL(i)=RecalclogLSprtParas(n,tmax,I,RL,tEs(:,k),tI,tRorD,tRL,tRLR,preB,tDm,prevK,IPm(:,k),r1,p1,h0,d,d0,beta,alpha,epsilon,delta,typ,tLm,i1,i2,i3,j1,j2,j3);
end

D=-2*mean(LL);