function [M,HPDI]=CalcWthnHHriskIncr(zthin,p,K0)
pthin=p(zthin,:);
K0thin=K0(zthin);
WHHRI=pthin(:,4)./(pthin(:,1).*K0thin);
figure;
[~,HPDI]=PlotPstrDistn(WHHRI,'Within-HH risk increase',50);
M=median(WHHRI);