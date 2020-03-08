function [mode_p,HPDI,mode_p1,HPDI1,pcorr,mode_sptl,HPDI_sptl,mode_bckgrnd,HPDI_bckgrnd,mode_d_half,HPDI_d_half,mode_WHHRI,HPDI_WHHRI,mean_IP,HPDI_IP]=ProcessOutput(str,burnin1,thin,doPlots,savePlots,sprtParas)
close all

% Load MCMC output
load(str)

% Discard burnins
if nargin==1
    z=burnin+1:niters; % if no burnin supplied, use saved default (round(niters/10))
    thin=1;
    doPlots=false;
    savePlots=false;
    sprtParas='';
else
    z=burnin1+1:niters; % if new burnin supplied, use it
end

%% TRACE PLOTS AND POSTERIOR DISTRIBUTIONS OF LOG-LIKELIHOOD AND PARAMETERS
nbins=50;
scrnsz=get(0,'ScreenSize');
zthin=z(1:thin:end);
figure;
% Uncomment line below to plot output with fixed parameters excluded and 1
% realisation of missing onset times
% [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(zthin,LL,p,nu,pname,prior_mean,p1,a,b,n,tmax,I,RpreD,DpreR,tI,tR,tD,tRLm,tRLRm,nbins,scrnsz);
% Plot output with multiple realisations of missing onset times
[mode_pu,HPDIu,mode_p1u,HPDI1u]=PlotOutput2(zthin,LL,p,nu,pname,prior_mean,p1,a,b,n,tmax,I,RpreD,DpreR,OR,NONR,RNO,ONR,tI,tR,tD,tRLm,tRLRm,tIsNONR,tIsRNO,tRsNONR,tRsONR,niters,nbins,scrnsz);
saveas2(gcf,['PSTR_DISTNS_' rslts],savePlots)
figure;
PlotTrace(zthin,p,nu,pname,p1,mode_pu,HPDIu,mode_p1u,HPDI1u,scrnsz)

figure;
% [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(zthin,LL,p,np,pname,prior_mean,p1,a,b,n,tmax,I,RpreD,DpreR,tI,tR,tD,tRLm,tRLRm,nbins,scrnsz);
[mode_p,HPDI,mode_p1,HPDI1]=PlotOutput2(zthin,LL,p,np,pname,prior_mean,p1,a,b,n,tmax,I,RpreD,DpreR,OR,NONR,RNO,ONR,tI,tR,tD,tRLm,tRLRm,tIsNONR,tIsRNO,tRsNONR,tRsONR,niters,nbins,scrnsz);
figure;
PlotTrace(zthin,p,np,pname,p1,mode_p,HPDI,mode_p1,HPDI1,scrnsz)

%% AUTOCORRELATION AND PARAMETER CORRELATION
% Autocorrelation plots for each parameter
figure; set(gcf, 'Position', [0 70 round(scrnsz(3)/2) round(scrnsz(3)/2.5)]);
for i=1:nu
    j=u(i);
    subplot(ceil((nu+1)/2),2,i)
    acf(p(z,j),min(200,numel(z)-1));
    title(['\' pname{j}])
end
subplot(ceil((nu+1)/2),2,nu+1)
% figure;
acf(p1(z),min(200,numel(z)-1));
title('$$p$$','Interpreter','latex')
saveas2(gcf,['ACFs_' rslts],savePlots)
saveas2(gcf,['ACFs_' rslts '.eps'],savePlots,'epsc')

% Matrix of parameter pair plots
figure; set(gcf, 'Position', [0 70 round(scrnsz(3)/2) round(scrnsz(3)/2.5)]);
[~,ax,~,~,~]=plotmatrix2([p(z,u),p1(z)]);
for i=1:nu
    j=u(i);
    xlabel(ax(end,i),['\' pname{j}],'Fontsize',13)
    ylabel(ax(i,1),['\' pname{j} '    '],'Fontsize',13,'rot',0)
end
xlabel(ax(end,nu+1),'$$p$$','Fontsize',13,'Interpreter','latex')
ylabel(ax(nu+1,1),'$$p$$    ','Fontsize',13,'Interpreter','latex','rot',0)
saveas2(gcf,['ParamCrrltn_' rslts],savePlots)
saveas2(gcf,['ParamCrrltn_' rslts '.eps'],savePlots,'epsc')

% Parameter correlation coefficients
pcorr=corrcoef([p(z,u) p1(z)]);

%% ESTIMATED TRANSMISSION KERNEL
% Plot estimated transmission kernel
if strcmp(sprtParas,'SprtParas')
    [mode_Ke,mode_K0,mode_rate,HPDI_rate]=PlotKnlSprtParas(zthin,p,K0,mode_p,HPDI,d,typ,n,i1,i2,i3,j1,j2,j3);
else
    [mode_Ke,mode_K0,mode_rate,HPDI_rate]=PlotKnl(zthin,p,K0,mode_p,HPDI,d,typ,n);
end
saveas2(gcf,['SPTL_KRNL_' rslts],savePlots)
saveas2(gcf,['SPTL_KRNL_' rslts '.eps'],savePlots,'epsc')

% Calculate transmission rate from VL cases and background transmission rate
mode_sptl=mode_rate(1)*1e4; % cases/10,000 people/mnth
HPDI_sptl=HPDI_rate*1e4; % cases/10,000 people/mnth
mode_bckgrnd=mode_p(3)*1e4; % cases/10,000 people/mnth
HPDI_bckgrnd=HPDI(3,:)*1e4; % cases/10,000 people/mnth

% Calculate half-risk distance analytically
if mode_p(2)~=0
    d_half=HalfRskDstnce(p(zthin,:),K0(zthin),typ);
    figure; 
    [mode_d_half,HPDI_d_half]=PlotPstrDistn(d_half,'d_{1/2}',200);
else
    mode_d_half=NaN;
    HPDI_d_half=NaN(1,2);
end

% Calculate within-HH risk increase (WHHRI)
if mode_p(4)==0
    mode_WHHRI=NaN;
    HPDI_WHHRI=NaN(1,2);
else
    [mode_WHHRI,HPDI_WHHRI]=CalcWthnHHriskIncr(zthin,p,K0);
end

%% INFECTION TIMES
if doPlots
    mode_tE=zeros(nI,1);
    nplot=20;
    tEc=NaN(nI,niters);
    for j=1:niters
        tEc(ismember(I,pick(:,j)),j)=tEs(ismember(I,pick(:,j)),j);
    end    
    % Cases with both onset and treatment times
    for i=1:nplot
        j=find(I==OR(i));
        figure;
        mode_tE(j)=PlotInfctnTimePstrDistn(tEc(j,z),tI(OR(i)),r1,p10,j);
        saveas2(gcf,['E' num2str(j) rslts],savePlots)
        saveas2(gcf,['E' num2str(j) rslts '.eps'],savePlots,'epsc')
    end
    %%
    % Cases without onset or treatment times
    for i=1:nplot
        j=find(I==NONR(i));
        figure;
        histogram(tEs(j,z),'Normalization','pdf','BinMethod','integers'); hold on
        histogram(tIsNONR(i,z),'Normalization','pdf','BinMethod','integers');
        histogram(tRsNONR(i,z),'Normalization','pdf','BinMethod','integers')
        set(gca,'FontSize',16);
        xlabel('t (months)','FontSize',16)
        ylabel('Density','FontSize',16)
        h1=legend(['$$E_{' num2str(j) '}$$'],['$$I_{' num2str(j) '}$$'],['$$R_{' num2str(j) '}$$']);
        set(h1,'Interpreter','latex')
        saveas2(gcf,['EIR' num2str(j) rslts],savePlots)
        if savePlots
            saveaspdf(gcf,['EIR' num2str(j) rslts])
        end
    end
    % Cases without onset times
    for i=1:nRNO
        j=find(I==RNO(i));
        figure;
        histogram(tEs(j,z),'Normalization','pdf','BinMethod','integers'); hold on
        histogram(tIsRNO(i,z),'Normalization','pdf','BinMethod','integers');
        set(gca,'FontSize',16);
        xlabel('t (months)','FontSize',16)
        ylabel('Density','FontSize',16)
        h2=legend(['$$E_{' num2str(j) '}$$'],['$$I_{' num2str(j) '}$$']);
        set(h2,'Interpreter','latex')
        saveas2(gcf,['EI' num2str(j) rslts],savePlots)
        if savePlots
            saveaspdf(gcf,['EI' num2str(j) rslts])
        end
    end
    %%
    mode_tR_ONR=zeros(nONR,1);
    % Cases without treatment times
    for i=1:nONR
        j=find(I==ONR(i));
        figure; 
        [mode_tE(j),hE]=PlotInfctnTimePstrDistn(tEs(j,z),tI(ONR(i)),r1,p10,ONR(i));
        hold on        
        [mode_tR_ONR(i),hR]=PlotRcvryTimePstrDistn(tRsONR(i,z),tI(ONR(i)),r0,p0,ONR(i));
        set(gca,'FontSize',16);
        xlabel('t (months)','Interpreter','tex','FontSize',16)
        ylabel('Density','FontSize',16)
        xlim([min(hE.BinEdges) max(hR.BinEdges)])
        h3=legend([hE hR],['$$E_{' num2str(j) '}$$'],['$$R_{' num2str(j) '}$$']);
        set(h3,'Interpreter','latex')
        saveas2(gcf,['ER' num2str(j) rslts],savePlots)
        if savePlots
            saveaspdf(gcf,['ER' num2str(j) rslts])
        end
    end
end
%% INCUBATION PERIODS
% Mean IP based on IP distn parameter
mean_IP=mean(r1*(1-p1(zthin))./p1(zthin))+1;
figure;
[mode_IP,HPDI_IP]=PlotPstrDistn(r1*(1-p1(zthin))./p1(zthin)+1,'IPD',50);

% Calculate vector of mean IPs for MCMC samples
mean_IPs=mean(IPs,1);
% Plot auto-correlation fn for mean incubation period
figure;
acf(mean_IPs(z)',min(200,numel(z)-1));
title('mean IP')

% Plot correlation between mean incubation period and p1
figure; plot(p1(z),mean_IPs(z),'.')
xlabel('p'); ylabel('mean IP')