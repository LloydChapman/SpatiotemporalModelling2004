function VLStmMCMC(db,r1,a,b,u,beta0,alpha0,epsilon0,delta0,h0,typ,inclLST,niters,plotOutpt,rslts)

if nargin==0
    db='data_final.mat'; % database
    r1=2; % r parameter of NegBin(r1,p1) distribution for incubation period
    % Shape parameters for beta prior for p1 parameter of NegBin(r1,p1) distn
    a=20; 
    b=38;
    u=1:3; % transmission parameters to update
    beta0=0.2; % spatial transmission rate constant
    alpha0=50; % distance scale factor for spatial kernel;
    epsilon0=1e-3; % background transmission rate
    delta0=0; % additional within-HH transmission rate
    h0=0.03; % relative infectiousness of pre-symptomatics
    typ='Cauchy'; % type of transmission kernel
    inclLST=false; % include LST data: false=don't include, true=include
    niters=2e5; % no. of MCMC steps
    plotOutpt=true; % plot output or not (if =false)
    rslts='MCMC_NBIP.mat'; % name of results file for output
end

%% LOAD DATA
load(db)

% Number of individuals
n=size(data,1);

% Set mean incubation period and treatment duration
mu=5; % months
% durRX=1; % month

% Set start year as 1998 and start month as Jan so infection times of those
% that had KA onset in Jan 1999 can be included
START_YR=1998;
START_MO=1;

%% MAKE VECTORS OF EVENT TIMES
% Find time of last event and make time vector
tmax=max(max(table2array(data(:,{'FEV_ONS_MOS','DX_MOS','RX_MOS','INT02_MOS','INT03_MOS','INT04_MOS','RELAPSE_MOS','RELP_RX_MOS','BIRTH_MOS','DEATH_MOS','LSTR02_MOS','LSTR03_MOS','LSTR04_MOS'}))));

% Make vectors of onset, recovery, birth and death times
tI=data.FEV_ONS_MOS;
tR=data.RX_MOS;
tRL=data.RELAPSE_MOS;
tRLR=data.RELP_RX_MOS;
tB=data.BIRTH_MOS;
tD=data.DEATH_MOS;
% Make index vectors for KA cases with and without onset times (KO and KNO), 
% treated indvdls, cases missing treatment times (NR), 2002 LST positives, 
% relapsers, treated relapsers, and births and deaths
OR=find(~isnan(tI)&(~isnan(tR)|~isnan(tD))); % KA: onset and treatment/death time
NONR=find(data.KA&isnan(tI)&isnan(tR)); % KA: no onset or treatment time
ONR=find(~isnan(tI)&isnan(tR)&isnan(tD)); % KA: onset but no treatment time & didn't die from KA
RNO=find(isnan(tI)&~isnan(tR)); % KA: treatment but no onset time
NO=[NONR;RNO];
I=sort([OR;NONR;ONR;RNO]);
NI=(isnan(data.PREVKA)&data.KA==0);
if inclLST
    L02=find(data.LST02==1&NI&isnan(tB)); % LST+ in 2002, excluding LST+ves who have had KA, 1 individual who has KA and is LST+, and 1 individual who gets KA later
    L02B=find(data.LST02==1&NI&~isnan(tB)); % LST+ in 2002, born after START time
    L03C=find(data.LST03==1&data.LST02==0&NI); % LST conversion between 2002 and 2003
    L03=find(data.LST03==1&isnan(data.LST02)&NI&isnan(tB)); % LST+ in 2003, no LST 2002
    L03B=find(data.LST03==1&isnan(data.LST02)&NI&~isnan(tB)); % LST+ in 2003, no LST in 2002, born after START time
    L04C=find(data.LST04==1&data.LST02==0&data.LST03==0&NI); % LST conversion between 2003 and 2004
    % N.B. No pstve LST04 asymptomatic individuals with missing LST02 or LST03
    L=sort([L02;L02B;L03C;L03;L03B;L04C]);
end
RL=find(~isnan(tRL));
RLR=find(~isnan(tRLR));
B=find(tB>0); % exclude people estimated to have been born before START time
D=find(~isnan(tD));
DpreR=find(~isnan(tI)&isnan(tR)&~isnan(tD));
RpreD=setdiff(I,DpreR);
% Number of KA cases
nI=numel(I);
% Number of cases with both onset and treatment times
nOR=numel(OR);
% Number of KA cases w/o onset or treatment times
nNONR=numel(NONR);
% Number of KA cases w/ onset time but no treatment time
nONR=numel(ONR);
% Number of KA cases w/ treatment time but no onset time
nRNO=numel(RNO);
% Number of cases with no onset time
nNO=numel(NO);

%% DRAW INITIAL MISSING ONSET AND TREATMENT TIMES
% Create vectors of lower and upper bounds for onset month
tIlb=NaN(n,1);
tIlb(I)=(data.FEV_ONS_YR(I)-START_YR)*12-START_MO+1+1;
tIub=NaN(n,1);
tIub(I)=(data.FEV_ONS_YR(I)-START_YR)*12-START_MO+1+12;
% Fit negative binomial distribution to onset-to-treatment (OT) times
[r0,p0]=FitOTdistn(tI,tR);
% For individuals with missing onset, diagnosis and treatment times draw
% onset time at random from onset year
for i=1:nNONR
    while isnan(tI(NONR(i))) || tI(NONR(i))>=tD(NONR(i))
        tI(NONR(i))=(data.FEV_ONS_YR(NONR(i))-START_YR)*12+randi(12)-START_MO+1;
    end
end
% Draw treatment times for these individuals using OT distn
for i=1:nNONR
    while isnan(tR(NONR(i))) || tR(NONR(i))>min(tD(NONR(i)),tmax)
        tR(NONR(i))=tI(NONR(i))+nbinrnd(r0,p0)+1;
    end
end
% For individuals with onset time but no treatment time draw treatment time
% using OT distn
for i=1:numel(ONR)
    while isnan(tR(ONR(i))) || tR(ONR(i))>min(tD(ONR(i)),tmax)
        tR(ONR(i))=tI(ONR(i))+nbinrnd(r0,p0)+1;
    end
end
% For individuals with treatment time but no onset time draw onset time at
% random from onset year
for i=1:numel(RNO)
    while isnan(tI(RNO(i))) || tI(RNO(i))>=tR(RNO(i))
        tI(RNO(i))=(data.FEV_ONS_YR(RNO(i))-START_YR)*12+randi(12)-START_MO+1;
    end
end

% Make vector of recovery times for treated individuals and death times for
% untreated individuals
tRorD=tR;
tRorD(DpreR)=tD(DpreR);

%% DRAW INITIAL INFECTION TIMES
% Draw infection times from negative binomial distribution with parameters
% r1 and p1
p10=r1/(mu-1+r1);
tE=NaN(n,1); % initialise infection time vector
tE(I)=0; % set infection times before start date initially
% Make vector of birth times for KA cases born during study
tBI=NaN(n,1);
tBI(setdiff(I,B))=0; % for those not born during study set dummy birth time of 0
tBI(intersect(I,B))=tB(intersect(I,B));
for i=1:nI
    % loop until infection time is after month 1 or birth month if
    % individual was born during study
    while tE(I(i))<max(1,tBI(I(i))+1) % can only be infected at least 1 month after birth
        tE(I(i))=tI(I(i))-(nbinrnd(r1,p10)+1);
    end
end

if inclLST
    if ~exist('LSTCnvsnTimes.mat','file')
        %% DRAW LST CONVERSION TIMES
        % Calculate probability that LST+ individuals at first skin test were LST+ at
        % START time using estimated conversion rate and asymptote from catalytic model
        [pars,~]=FitCatModLST(data);
        l=pars(1);
        c=pars(2);
        % Time first skin test read (in months)
        tLR=[data.LSTR02_MOS(L02);data.LSTR03_MOS(L03)];
        % Ages of LST positive individuals when first skin test read
        a1=[data.AGE02(L02);data.AGE03(L03)];
        % Ages of LST positive individuals at START time
        a0=a1-tLR/12;
        
        % Probability of being LST+ at START time based on age
        P_Lpos0=PrevVA(a0,l,c);
        % Probability of being LST+ in 2002 based on age
        P_Lpos1=PrevVA(a1,l,c);
        % Probability of being LST+ at START time given LST+ in 2002
        P_Lpos0pos1=P_Lpos0./P_Lpos1;

        % Use probabilities to set whether individuals were LST+ or LST- at
        % START time
        L0203=[L02;L03];
        Lpos0=L0203(P_Lpos0pos1>rand(numel(L0203),1));
        Lneg0=setdiff(L0203,Lpos0);
        % Make vector of LST conversion dates
        tL=NaN(n,1);
        tL(Lpos0)=1;
        % Draw LST conversion month from exponential distribution for those LST-
        % at START time
        for i=1:numel(Lneg0)
            % resample until LST conversion is before (or in same month as) first test
            while isnan(tL(Lneg0(i))) || tL(Lneg0(i))>tLR(L0203==Lneg0(i))
                tL(Lneg0(i))=ceil(exprnd(12/l));
            end
        end
        % Draw LST conversion date for those born after START time who are LST+ in 2002
        for i=1:numel(L02B)
            % resample until LST conversion is before (or in same month as) 2002 test
            while isnan(tL(L02B(i))) || tL(L02B(i))>data.LSTR02_MOS(L02B(i)) 
                tL(L02B(i))=tB(L02B(i))+ceil(exprnd(12/l));
            end
        end
        % Draw LST conversion date for those born after START time who are first measured LST+ in 2003
        for i=1:numel(L03B)
            % resample until LST conversion is before (or in same month as) 2003 test
            while isnan(tL(L03B(i))) || tL(L03B(i))>data.LSTR03_MOS(L03B(i)) 
                tL(L03B(i))=tB(L03B(i))+ceil(exprnd(12/l));
            end
        end
        % Draw time of LST conversion after 1st LST test for those who became LST+ between 1st and 2nd test
        for i=1:numel(L03C)
            % resample until LST conversion is before (or in same month as) 2003 test
            while isnan(tL(L03C(i))) || tL(L03C(i))>data.LSTR03_MOS(L03C(i)) 
                tL(L03C(i))=data.LSTR02_MOS(L03C(i))+ceil(exprnd(12/l));
            end
        end
        % Draw time of LST conversion after 2nd LST test for those who became LST+ between 2nd and 3rd test
        for i=1:numel(L04C)
            % resample until LST conversion is before (or in same month as) 2004 test
            while isnan(tL(L04C(i))) || tL(L04C(i))>data.LSTR04_MOS(L04C(i))
                tL(L04C(i))=data.LSTR03_MOS(L04C(i))+ceil(exprnd(12/l));
            end
        end
        % Set any LST conversion times before START time (for individuals with tB<1) to 1
        tL(tL<1)=1;
        save('LSTCnvsnTimes','tL');
    else
        load('LSTCnvsnTimes.mat');
    end
end

%% MAKE STATUS MATRICES
% Make logical matrices for infection, relapse, relapse recovery, birth and death times
% Infection
tEm=false(n,tmax);
tEm((tE(I)-1)*n+I)=1;
if inclLST
    % LST conversion
    tLm=false(n,tmax);
    tLm((tL(L)-1)*n+L)=1;
end
% Relapse
tRLm=false(n,tmax);
tRLm((tRL(RL)-1)*n+RL)=1;
% Relapse recovery
tRLRm=false(n,tmax);
tRLRm((tRLR(RLR)-1)*n+RLR)=1;
% Birth
tBm=false(n,tmax);
tBm((tB(B)-1)*n+B)=1;
% Death
tDm=false(n,tmax);
tDm((tD(D)-1)*n+D)=1;
% Pre-birth
preB=false(n,tmax);
preB(B,:)=1-cumsum(tBm(B,:),2);

% Index vector for previous KA
prevK=(data.PREVKA==1);

% Make susceptibility status matrix
if inclLST
    S=1-preB-max(max(cumsum(tEm,2),cumsum(tDm,2)),cumsum(tLm,2)); % remove LST+ individuals from susceptibles
else
    S=1-preB-max(cumsum(tEm,2),cumsum(tDm,2)); % don't remove LST+ individuals
end
S(prevK,:)=0; % remove previous KA cases from susceptibles

%% SET UP MCMC
burnin=round(niters/10); % burn-in
nu=numel(u); % number of parameters to update
pname={'beta','alpha','epsilon','delta'}; % parameter labels
np=numel(pname); % number of possible parameters in model

% Calculate initial asymptomatic infection periods
IPs=zeros(nI,niters); % matrix of infection periods
IPold=NaN(n,1);
IPold(I)=tI(I)-tE(I);

% Infectiousness over time
% Assume constant relative infectiousness h0 from infection to onset and 1 
% from onset to recovery
h=zeros(n,tmax); % size of infectiousness matrix
for j=1:nI
    h(I(j),tE(I(j)):tI(I(j))-1)=h0; % infectiousness from month of infection up to month before symptom onset
    h(I(j),tI(I(j)):tRorD(I(j))-1)=1; % infectiousness up to month before month of treatment or death
end
% Assume individuals return to full infectiousness on relapse until
% re-treatment
for j=1:numel(RL)
    h(RL(j),tRL(RL(j)):min(tRLR(RL(j))-1,tmax))=1;
end

% Calculate distances between HHs
d=CalcDists(data);

% Calculate infectious pressure over time
[K,K0old]=Knl(d,alpha0,typ,n); % spatial kernel
d0=double(d==0); % indicator matrix of individuals living in the same HHs
d0(1:n+1:end)=0; % set diagonal (same individual) entries to 0
rate=beta0*K+delta0*d0; % calculate transmission rate, adding within-HH contribution
lambda=rate(:,I)*h(I,:)+epsilon0; % infectious pressure

% PRIORS
prior_mean=[1,50,1,1]; % prior distribution means

% PROPOSALS
ppmean=zeros(1,np); % vector for storing means for proposal distributions for transmission parameters
ppvar0=diag([beta0,alpha0,epsilon0,delta0].^2); % initial variances for proposal distribution for block update
ppvar=zeros(np); % vector for storing variances for proposal distributions for transmission parameters
nEmoves=round(nOR/5); % number of cases with onset and recovery times to propose new infection times for in each step 
pick=zeros(nEmoves+nNONR+nONR+nRNO,niters);
for i=1:niters 
pick(:,i)=[OR(randperm(nOR,nEmoves));NONR;ONR;RNO]; % indices of individuals to propose new infection times for in each step
% N.B. randi rather than randperm to avoid possibility of changing same infection time twice in one step
end
ERvar=4; % variance for infection period moves

% INITIALISATION
% Fitting 4 (or 5) transmission parameters: beta, alpha, epsilon (and delta) and p1
pold=[beta0,alpha0,epsilon0,delta0]; % set old parameter values
plb=zeros(1,np); % lower limits for parameters
pub=Inf(1,np); % upper limits for parameters
p1new=p10; % initial p1 value
LLold=logL(S,lambda,tEm,IPold(I),r1,p10); % initial log-likelihood

% Matrices and vectors for saving parameter and log-likelihood values
p=zeros(niters,np); % matrix for saving parameter values
p1=zeros(niters,1); % vector for saving p1 values
K0=zeros(niters,1); % vector for saving spatial kernel normalisation constants
LL=zeros(niters,1); % vector for saving log-likelihood values
terms=zeros(niters,3); % saving individual log-likelihood terms
tEs=zeros(nI,niters);
tIsNONR=zeros(nNONR,niters);
tRsNONR=zeros(nNONR,niters);
tIsRNO=zeros(nRNO,niters);
tRsONR=zeros(nONR,niters);

% Initialise new status and infectiousness matrices
pnew=pold;
tEnew=tE;
tInew=tI;
tRnew=tR;
tEmnew=tEm;
Snew=S; 
IPnew=IPold;
hnew=h;

% Initialise acceptance and rejection counts and rates
acc_p=0; % initial acceptance count for transmission parameter updates
rej_p=0; % initial rejection count for transmission parameter updates
acc_rate_p=0; % initial acceptance rate for transmission parameter updates
acc_E=0;
rej_E=0;
acc_rate_E=0;
acc_I=0;
rej_I=0;
acc_rate_I=0;
acc_ERmove=0;
rej_ERmove=0;
acc_rate_ERmove=0;
acc_R=0;
rej_R=0;
acc_rate_R=0;

nbins=50;
scrnsz=get(0,'ScreenSize');

for k=1:niters
    %% UPDATE TRANSMISSION PARAMETERS USING ADAPTIVE RANDOM WALK METROPOLIS-HASTINGS    
    if k<=100 || rand<0.05
        pnew(u)=mvnrnd(pold(u),0.1^2*ppvar0(u,u)/nu);
        while any(pnew(u)<=plb(u)) || any(pnew(u)>=pub(u)) % resample if outside parameter range
            pnew(u)=mvnrnd(pold(u),0.1^2*ppvar0(u,u)/nu);
        end
    else
        pnew(u)=mvnrnd(pold(u),2.38^2*ppvar(u,u)/nu);
        while any(pnew(u)<=plb(u)) || any(pnew(u)>=pub(u)) % resample if outside parameter range
            pnew(u)=mvnrnd(pold(u),2.38^2*ppvar(u,u)/nu);
        end
    end
    
    % Calculate ratio of prior probabilities for new and old value
    q=log(prod(exppdf(pnew(u),prior_mean(u))))-log(prod(exppdf(pold(u),prior_mean(u))));
    [K,K0new]=Knl(d,pnew(2),typ,n); % update distance kernel
    rate_new=pnew(1)*K+pnew(4)*d0; % update rate
    lambda_new=rate_new(:,I)*h(I,:)+pnew(3); % update whole infectious pressure
    
    [LLnew,LL1,LL2,LL3]=logL(S,lambda_new,tEm,IPold(I),r1,p1new); % calculate log-likelihood
    
    log_ap=LLnew-LLold+q; % calculate Metropolis-Hastings acceptance probability
    
    if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
        pold=pnew; % keep updated parameter values
        lambda=lambda_new; rate=rate_new; K0old=K0new; LLold=LLnew; % keep updated info
        acc_p=acc_p+1; % add to acceptance
    else
        lambda_new=lambda; % revert lambda_new to old value as only part of lambda_new is updated in infection time updates below
        rej_p=rej_p+1; % add to rejection
    end
    
    %% GIBBS UPDATE SUCCESS PROBABILITY PARAMETER FOR ASYMPTOMATIC INFECTION PERIOD DISTRIBUTION
    p1new=betarnd(a+r1*nI,b+sum(IPold(I))-nI);

    %% UPDATE INFECTION TIMES    
    for i=1:size(pick,1)
        j=pick(i,k); % get index of infection time to update
        tEnew(j)=0; % set infection time as before start date initially
        while tEnew(j)<max(1,tBI(j)+1) % can only be infected at least 1 month after birth
            tEnew(j)=tI(j)-(nbinrnd(r1,p1new)+1); % sample from negative binomial distribution
        end
        tEj=tE(j);
        tEjnew=tEnew(j);
        IPnew(j)=tI(j)-tEjnew; % new incubation period
        q=log(nbinpdf(IPold(j)-1,r1,p1new))-log(nbinpdf(IPnew(j)-1,r1,p1new)); % proposal ratio for infection time update
       
        tEmnew(j,tEj)=0; % remove old infection time
        tEmnew(j,tEjnew)=1; % add new infection time
        Snew(j,tEjnew:tEj-1)=0; % remove old susceptible times if new infection time is earlier (from month indvdl is infctd up to but not incl. month before old infection time)
        Snew(j,tEj:tEjnew-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)
        
        hnew(j,tEj:tEjnew-1)=0; % remove old infectiousness if new infection time is later
        hnew(j,tEjnew:tEj-1)=h0; % add infectiousness if new infection time is earlier
        erlrE=min(tEj,tEjnew); % index of column for earlier infctn time between old and new infection time
        ltrE=max(tEj,tEjnew); % index of column for later infctn time between old and new infection time
        lambda_new(:,erlrE:ltrE-1)=rate(:,I)*hnew(I,erlrE:ltrE-1)+pold(3); % update infectious pressure
        
        [LLnew,LL1,LL2,LL3]=logL(Snew,lambda_new,tEmnew,IPnew(I),r1,p1new); % calculate log-likelihood
        
        log_ap=LLnew-LLold+q; % calculate Metropolis-Hastings acceptance probability
        
        if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
            IPold(j)=IPnew(j); tE(j)=tEnew(j); S(j,:)=Snew(j,:); tEm(j,:)=tEmnew(j,:);
            h(j,:)=hnew(j,:); lambda(:,erlrE:ltrE-1)=lambda_new(:,erlrE:ltrE-1); LLold=LLnew; % keep updated info, overwrite old log-likelihood for step k with new log-likelihood
            acc_E=acc_E+1;
        else
            IPnew(j)=IPold(j); tEnew(j)=tE(j); Snew(j,:)=S(j,:); tEmnew(j,:)=tEm(j,:);
            hnew(j,:)=h(j,:); lambda_new(:,erlrE:ltrE-1)=lambda(:,erlrE:ltrE-1); % keep old values, don't change log-likelihood
            rej_E=rej_E+1;
        end
    end
    
    %% UPDATE MISSING ONSET TIMES
    for i=1:nNO
        j=NO(i);
        tInew(j)=0;     
        while tInew(j)<=tE(j) || tInew(j)<tIlb(j) || tInew(j)>tIub(j) || tInew(j)>=tD(j)
            tInew(j)=tR(j)-(nbinrnd(r0,p0)+1);
        end
        tIj=tI(j);
        tIjnew=tInew(j);
        IPnew(j)=tIjnew-tE(j); % new incubation period
        hnew(j,tIj:tIjnew-1)=h0; % reduce infectiousness up to new onset time if it is later
        hnew(j,tIjnew:tIj-1)=1; % increase infectiousness from new onset time if it is earlier
        erlrI=min(tIj,tIjnew); % index of column for earlier onset time between old and new onset time
        ltrI=max(tIj,tIjnew); % index of column for later onset time between old and new onset time
        lambda_new(:,erlrI:ltrI-1)=rate(:,I)*hnew(I,erlrI:ltrI-1)+pold(3); % update infectious pressure
        
        [LLnew,LL1,LL2,LL3]=logL(S,lambda_new,tEm,IPnew(I),r1,p1new); % calculate log-likelihood
        
        log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability

        if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
            IPold(j)=IPnew(j); tI(j)=tInew(j); h(j,:)=hnew(j,:); lambda(:,erlrI:ltrI-1)=lambda_new(:,erlrI:ltrI-1); LLold=LLnew; % keep updated info, overwrite old log-likelihood for step k with new log-likelihood
            acc_I=acc_I+1;
        else
            IPnew(j)=IPold(j); tInew(j)=tI(j); hnew(j,:)=h(j,:); lambda_new(:,erlrI:ltrI-1)=lambda(:,erlrI:ltrI-1); % keep old values, don't change log-likelihood
            rej_I=rej_I+1;
        end
    end

    %% MOVE WHOLE INFECTION-TO-TREATMENT PERIOD
    for i=1:nNONR
        j=NONR(i);
        tmp=0;
        % if onset and death in same year and maximum symptom period 
        % currently chosen, or infection is at start of study/birth time
        % and onset/recovery time are as late as they can be, do nothing as
        % infection-to-treatment period cannot be moved
        if ~(tI(j)==tIlb(j) && tR(j)==tD(j)) && ~(tE(j)==max(1,tBI(j)+1) && (tI(j)==min(tIub(j),tD(j)-1) || tR(j)==tD(j))) % only propose infection-treatment block move if a move is possible
            % resample until different onset time (in onset year and before
            % death month) and treatment time (before death month) are chosen;
            % allow treatment month to be same as death month as many cases
            % were said to have died during treatment
            while tmp==0 || tInew(j)<tIlb(j) || tInew(j)>tIub(j) || tInew(j)>=tD(j) || tEnew(j)<max(1,tBI(j)+1) || tRnew(j)>min(tD(j),tmax)
                tmp=round(sqrt(ERvar)*randn);
                tEnew(j)=tE(j)+tmp;
                tInew(j)=tI(j)+tmp;
                tRnew(j)=tR(j)+tmp;
            end
            tEj=tE(j);
            tEjnew=tEnew(j);
            tEmnew(j,tEj)=0; % remove old infection time
            tEmnew(j,tEjnew)=1; % add new infection time
            Snew(j,tEjnew:tEj-1)=0; % remove old susceptible times if new infection time is earlier (from month indvdl is infctd up to but not incl. month before old infection time)
            Snew(j,tEj:tEjnew-1)=1; % add new susceptible times if new infection time is later (up to but not incl. month indvdl is infctd)
            
            hnew(j,tEj:tEjnew-1)=0; % remove old infectiousness if new infection time is later
            hnew(j,tEjnew:tEj-1)=h0; % add infectiousness if new infection time is earlier
            erlrE=min(tEj,tEjnew); % index of column for earlier infctn time between old and new infection time
            ltrE=max(tEj,tEjnew); % index of column for later infctn time between old and new infection time
            
            tIj=tI(j);
            tIjnew=tInew(j);
            hnew(j,tIj:tIjnew-1)=h0; % reduce infectiousness up to new onset time if it is later
            hnew(j,tIjnew:tIj-1)=1; % increase infectiousness from new onset time if it is earlier
            erlrI=min(tIj,tIjnew); % index of column for earlier onset time between old and new onset time
            ltrI=max(tIj,tIjnew); % index of column for later onset time between old and new onset time
            
            tRj=tR(j);
            tRjnew=tRnew(j);
            hnew(j,tRj:tRjnew-1)=1; % increase infectiousness up to new treatment time if it is later
            hnew(j,tRjnew:tRj-1)=0; % reduce infectiousness from new treatment time if it is later
            erlrR=min(tRj,tRjnew); % index of column for earlier treatment time between old and new treatment time
            ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
            
            lambda_new(:,[erlrE:ltrE-1,erlrI:ltrI-1,erlrR:ltrR-1])=rate(:,I)*hnew(I,[erlrE:ltrE-1,erlrI:ltrI-1,erlrR:ltrR-1])+pold(3); % update infectious pressure

            [LLnew,LL1,LL2,LL3]=logL(Snew,lambda_new,tEmnew,IPold(I),r1,p1new); % calculate log-likelihood
            
            log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance probability
            
            if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
                tE(j)=tEnew(j); S(j,:)=Snew(j,:); tEm(j,:)=tEmnew(j,:); tI(j)=tInew(j); tR(j)=tRnew(j);
                h(j,:)=hnew(j,:); lambda(:,[erlrE:ltrE-1,erlrI:ltrI-1,erlrR:ltrR-1])=lambda_new(:,[erlrE:ltrE-1,erlrI:ltrI-1,erlrR:ltrR-1]); LLold=LLnew; % keep updated info, overwrite old log-likelihood for step k with new log-likelihood
                acc_ERmove=acc_ERmove+1;
            else
                tEnew(j)=tE(j); Snew(j,:)=S(j,:); tEmnew(j,:)=tEm(j,:); tInew(j)=tI(j); tRnew(j)=tR(j);
                hnew(j,:)=h(j,:); lambda_new(:,[erlrE:ltrE-1,erlrI:ltrI-1,erlrR:ltrR-1])=lambda(:,[erlrE:ltrE-1,erlrI:ltrI-1,erlrR:ltrR-1]); % keep old values, don't change log-likelihood
                rej_ERmove=rej_ERmove+1;
            end
        end
    end
    
    %% UPDATE MISSING TREATMENT TIMES
    for i=1:nONR
        j=ONR(i);
        tRnew(j)=NaN;       
        while isnan(tRnew(j)) || tRnew(j)>min(tD(j),tmax)
            tRnew(j)=tI(j)+nbinrnd(r0,p0)+1;
        end
        tRj=tR(j);
        tRjnew=tRnew(j);
        hnew(j,tRj:tRjnew-1)=1; % increase infectiousness up to new treatment time if it is later
        hnew(j,tRjnew:tRj-1)=0; % reduce infectiousness from new treatment time if it is later
        erlrR=min(tRj,tRjnew); % index of column for earlier treatment time between old and new treatment time
        ltrR=max(tRj,tRjnew); % index of column for later treatment time between old and new treatment time
        lambda_new(:,erlrR:ltrR-1)=rate(:,I)*hnew(I,erlrR:ltrR-1)+pold(3); % update infectious pressure

        [LLnew,LL1,LL2,LL3]=logL(S,lambda_new,tEm,IPold(I),r1,p1new); % calculate log-likelihood
        
        log_ap=LLnew-LLold; % calculate Metropolis-Hastings acceptance value

        if log_ap > log(rand) % for acc_prob<=1, change if greater than rand
            tR(j)=tRnew(j); h(j,:)=hnew(j,:); lambda(:,erlrR:ltrR-1)=lambda_new(:,erlrR:ltrR-1); LLold=LLnew; % keep updated info, overwrite old log-likelihood for step k with new log-likelihood
            acc_R=acc_R+1;
        else
            tRnew(j)=tR(j); hnew(j,:)=h(j,:); lambda_new(:,erlrR:ltrR-1)=lambda(:,erlrR:ltrR-1); % keep old values, don't change log-likelihood
            rej_R=rej_R+1;
        end
    end
    
    %% SAVE PARAMETER VALUES AND PLOT PROGRESS
    % Save parameters, log-likelihood, and asymptomatic infection periods
    p(k,:)=pold;
    p1(k)=p1new;
    K0(k)=K0old;
    LL(k)=LLold;
    terms(k,:)=[LL1,LL2,LL3];   
    IPs(:,k)=IPold(I);
    tEs(:,k)=tE(I);
    tIsNONR(:,k)=tI(NONR);
    tRsNONR(:,k)=tR(NONR);
    tIsRNO(:,k)=tI(RNO);
    tRsONR(:,k)=tR(ONR);
    
    % Update empirical mean and covariance for proposal distribution
    [ppmean,ppvar]=updateMeanAndCov(ppmean,ppvar,pold,k);
    
    % Calculate acceptance rates for parameters and infection time moves
    acc_rate_p=acc_p/(acc_p+rej_p);
    acc_rate_E=acc_E/(acc_E+rej_E);
    acc_rate_I=acc_I/(acc_I+rej_I);
    acc_rate_ERmove=acc_ERmove/(acc_ERmove+rej_ERmove);
    acc_rate_R=acc_R/(acc_R+rej_R);
    
    % Print output
    if mod(k,1e3)==0 || k==niters % every 1000 iterations
        fprintf('Iteration %d done.\n', k); % display iteration number
        fprintf('Current likelihood=%6.4g\n', LL(k));
        fprintf('LL1=%6.4g\n', terms(k,1));
        fprintf('LL2=%6.4g\n', terms(k,2));
        fprintf('LL3=%6.4g\n', terms(k,3));
        
        % Plot output
        if plotOutpt && k>burnin % ignore burn-in period
            z=burnin+1:k; % iterations to plot
            figure(4);
            [mode_p,HPDI,mode_p1,HPDI1]=PlotOutput(z,LL,p,np,pname,prior_mean,p1,a,b,n,tmax,I,RpreD,DpreR,tI,tR,tD,tRLm,tRLRm,nbins,scrnsz);            
            figure(5);
            PlotTrace(z,p,np,pname,p1,mode_p,HPDI,mode_p1,HPDI1,scrnsz)
            drawnow
            for j=1:np
                fprintf(['mode ' pname{j} '=%6.4g\n'], mode_p(j));
            end
            fprintf('mode p1=%6.4g\n', mode_p1);
        end
    end
end

% Display acceptance rates
fprintf('Parameter block update acceptance rate is %5.3f%%.\n',100*acc_rate_p);
fprintf('Infection time move acceptance rate is %5.3f%%.\n',100*acc_rate_E);
fprintf('Onset time move acceptance rate is %5.3f%%.\n',100*acc_rate_I);
fprintf('Infection period move acceptance rate is %5.3f%%.\n',100*acc_rate_ERmove);
fprintf('Treatment time move acceptance rate is %5.3f%%.\n',100*acc_rate_R);

save(rslts)
