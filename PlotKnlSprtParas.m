function [mode_Ke,mode_K0,mode_rate,HPDIrate]=PlotKnlSprtParas(zthin,p,K0,mode_p,HPDI,d,typ,n,i1,i2,i3,j1,j2,j3)
dd=0:500;
nz=numel(zthin);
pthin=p(zthin,:);
K0thin=K0(zthin);
HPDIlb=repmat(HPDI(:,1)',nz,1);
HPDIub=repmat(HPDI(:,2)',nz,1);
idx=all(pthin>=HPDIlb & pthin<=HPDIub,2);
Ke=NaN(sum(idx),numel(dd));
for j=1:numel(dd)
    Ke(:,j)=KnlEst(dd(j),pthin(idx,2),typ);
end
rate=bsxfun(@times,Ke,pthin(idx,1).*K0thin(idx));
rate(:,1)=rate(:,1)+pthin(idx,4);
rateLB=min(rate,[],1);
rateUB=max(rate,[],1);
HPDIrate=[rateLB(1) rateUB(1)];
mode_Ke=KnlEst(dd,mode_p(2),typ);
mode_K0=NrmlstnConstSprtParas(d,mode_p(2),typ,n,i1,i2,i3,j1,j2,j3);
mode_rate=mode_p(1)*mode_K0*mode_Ke;
mode_rate(1)=mode_rate(1)+mode_p(4);
dd2=[dd,fliplr(dd)];
HPDReps=[HPDI(3,1)*ones(1,numel(dd)),HPDI(3,2)*ones(1,numel(dd))];
figure; fill(dd2,HPDReps,[0.9 0.9 0.9],'LineStyle','none')
hold on
plot([0 0],[rateLB(1) rateUB(1)],'Color',[0.7 0.7 0.7],'LineWidth',2)
HPDR=[rateLB,fliplr(rateUB)];
fill(dd2,HPDR,[0.7 0.7 0.7],'LineStyle','none')
h1=plot(dd([1 end]),[mode_p(3) mode_p(3)],'k--','LineWidth',1.5);
h2=plot(dd,mode_rate,'k','LineWidth',2);
set(gca,'FontSize',16);
xlabel('d (m)','FontSize',20);
ylabel('Transmission rate (mnth^{-1})')
if mode_p(4)==0
    legend([h2 h1],'\beta K(d)','\epsilon');
elseif mode_p(4)~=0
    legend([h2 h1],'\beta K(d)+\delta 1_{d=0}','\epsilon');
end
hold off
