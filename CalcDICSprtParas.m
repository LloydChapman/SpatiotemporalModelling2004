function DIC=CalcDICSprtParas(rslts1,burnin1,ss,str)

load(rslts1)

if nargin==1
    burnin1=round(niters/10);
    ss='mode';
    str='DIC2';
end

if ~exist('z','var')
    z=burnin1+1:niters;
end

if ~exist('tLm','var')
    tLm=[];
end

if exist('IPs','var')
    [D,LLM]=devSprtParas(p,p1,z,n,tmax,NONR,RNO,ONR,tIsNONR,tRsNONR,tIsRNO,tRsONR,I,RL,tEs,tI,tRorD,tRL,tRLR,preB,tDm,prevK,IPs,r1,h0,d,d0,typ,tLm,ss,i1,i2,i3,j1,j2,j3);
else
    [D,LLM]=devSprtParas(p,p1,z,n,tmax,NONR,RNO,ONR,tIsNONR,tRsNONR,tIsRNO,tRsONR,I,RL,tEs,tI,tRorD,tRL,tRLR,preB,tDm,prevK,IPm,r1,h0,d,d0,typ,tLm,ss,i1,i2,i3,j1,j2,j3);
end
MD=-2*mean(LL(z));
DIC=dic(D,MD);

save([str '_' rslts1],'D','LLM','MD','DIC')