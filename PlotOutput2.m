function [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput2(z,LL,p,np,pname,prior_mean,p1,a,b,n,tmax,I,RpreD,DpreR,OR,NONR,RNO,ONR,tI,tR,tD,tRLm,tRLRm,tIsNONR,tIsRNO,tRsNONR,tRsONR,niters,nbins,scrnsz,u)
set(gcf, 'Position', [0 70 round(scrnsz(3)/2) scrnsz(4) - 150]);
mode_p=zeros(1,np); % modal parameter values
HPDI=zeros(np,2); % highest posterior density intervals
subplot(3+floor(np/2), 2, [1 2])
plot(z,LL(z));
axis([z(1) z(end) min(LL(z)) max(LL(z))]);
xlabel('Iteration');
ylabel('Log likelihood');
for j=1:np    
    subplot(3+floor(np/2), 2, 2+j)
    [mode_p(j),HPDI(j,:)]=PlotPstrDistn(p(z,j),pname{j},nbins,'exp',prior_mean(j));
end
subplot(3+floor(np/2),2,np+3)
[mode_p1,HPDI1]=PlotPstrDistn(p1(z),'p',nbins,'beta',a,b);
subplot(3+floor(np/2),2,[5+2*floor(np/2) 6+2*floor(np/2)])

% Number of iterations to plot
npp=min(100,niters);
pp=1:round(niters/npp):niters;

% Calculate numbers of KA cases over time for different MCMC iterations
tIs=repmat(tI,1,npp);
tIs(NONR,:)=tIsNONR(:,pp);
tIs(RNO,:)=tIsRNO(:,pp);
tIm=false(n,tmax,npp);
tRs=repmat(tR,1,npp);
tRs(NONR,:)=tRsNONR(:,pp);
tRs(ONR,:)=tRsONR(:,pp);
tRorDm=false(n,tmax,npp);
for i=1:npp
tIm((i-1)*n*tmax+(tIs(I,i)-1)*n+I)=1;
tRorDm((i-1)*n*tmax+(tRs(RpreD,i)-1)*n+RpreD)=1;
tRorDm((i-1)*n*tmax+(tD(DpreR)-1)*n+DpreR)=1;
end
tRLm_rep=repmat(tRLm,1,1,npp);
tRLRm_rep=repmat(tRLRm,1,1,npp);
numI=sum(cumsum(tIm,2)-cumsum(tRorDm,2)+cumsum(tRLm_rep,2)-cumsum(tRLRm_rep,2));

% Calculate numbers of KA cases over time with observed onset and recovery
tIm1=false(n,tmax);
tIm1((tI(OR)-1)*n+OR)=1;
ORpreD=setdiff(OR,DpreR);
ODpreR=intersect(OR,DpreR);
tRorDm1=false(n,tmax);
tRorDm1((tR(ORpreD)-1)*n+ORpreD)=1;
tRorDm1((tD(ODpreR)-1)*n+ODpreR)=1;
numI1=sum(cumsum(tIm1,2)-cumsum(tRorDm1,2)+cumsum(tRLm,2)-cumsum(tRLRm,2));

% Plot imputed case numbers vs observed case numbers
jj=13:tmax;
plot(jj,squeeze(numI(1,jj,:)),'r'); hold on
plot(jj,numI1(jj),'k--','LineWidth',2); hold off
xlim([jj(1) Inf])
xlabel('t (months)');
ylabel('No. of VL cases');