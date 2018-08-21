function [pars,CI]=FitCatModLST(data)
%FITCATMODLST Fit variable asymptote catalytic model to LST age prevalence data

%% CALCULATE PREVALENCE OF LST POSITIVITY FOR EACH AGE
% Indicator for no current or previous KA
NI=(isnan(data.PREVKA)&data.KA==0);
% Make vector of ages at which individuals were tested with LST
a=3:max(data.AGE02(~isnan(data.LST02)&NI));
% Sum number tested at each age
na=zeros(numel(a),1);
for i=1:numel(a)
    na(i)=sum(data.AGE02==a(i)&~isnan(data.LST02)&NI);
end
% Sum number of LST positives at each age
naLpos=zeros(numel(a),1);
for i=1:numel(a)
    naLpos(i)=sum(data.AGE02==a(i)&data.LST02==1&NI);
end

% Remove ages at which no one was tested
iL=(na~=0);
a=a(iL);
naLpos=naLpos(iL);
na=na(iL);
% Calculate prevalence
prevLpos=naLpos./na;

% figure; plot(a,-log(prevLpos),'x');
% figure; bar(a,prevLpos)

%% CATALYTIC MODELS
data1(1,:)=a;
data1(2,:)=a+1;
data1(3,:)=na;
data1(4,:)=naLpos;

% VARIABLE ASYMPTOTE (VA) CATALYTIC MODEL
% Set initial guesses for lambda and c
lambda0=0.1;
c0=0.75;
[pars,CI]=mle(data1(:),'nloglf',@(params,data,cens,freq)NegLLVA(params,data,cens,freq),'start',[lambda0,c0],'lowerbound',[0,0],'upperbound',[Inf,1]);

% Plot prevalence of LST positivity
aa=linspace(0,max(data.AGE02),1000);
p=PrevVA(aa,pars(1),pars(2));
figure; plot(a,prevLpos,'x',aa,p,'r',aa([1,end]),[pars(2) pars(2)],'k--','LineWidth',1.5);
legend('data','model','c');
xlabel('Age (years)');
ylabel('Proportion LST+')
saveas(gcf,'LSTAgePrev')
saveas(gcf,'LSTAgePrev.eps','epsc')

save('LSTCatModelOutput','pars','L','a','prevLpos','aa','p')