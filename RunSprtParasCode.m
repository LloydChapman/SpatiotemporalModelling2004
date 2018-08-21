clear
close all

%% Set models to run
%% Inputs
db='data_final.mat'; % database
r=3; % negative binomial incubation period
mu=5; % guess for mean incubation period
% Shape parameters for beta prior for p1
a=22;
b=(mu-1)*(a-1)/r; % since mu-1=r*b/(a-1)

% Initial guesses for transmission parameters
beta0=0.2; % spatial transmission rate constant
alpha0=50; % distance scale factor for spatial kernel;  
epsilon0=1e-3; % background transmission rate
delta0=1e-3; % additional within-HH transmission rate
h0=0.03; % relative infectivity of pre-symptomatics
niters=4e5; % number of MCMC iterations
plotOutpt=false; % don't plot output

if r>1
    IPD='NBIP';
elseif r==1
    IPD='GIP';
end

% Set variables for different models
Mdls={'_BCKGRND_ONLY','_NONSPTL','','_DELTA','_EXP_KNL','_EXP_KNL_DELTA','_PRESX_NONINFCTS','_DELTA_PRESX_NONINFCTS','_EXP_KNL_PRESX_NONINFCTS','_EXP_KNL_DELTA_PRESX_NONINFCTS','_EXP_KNL_HGHR_PRESX_INFCTS','_EXP_KNL_DELTA_HGHR_PRESX_INFCTS','_EXP_KNL_HGHST_PRESX_INFCTS','_EXP_KNL_DELTA_HGHST_PRESX_INFCTS','_LST','_DELTA_LST','_EXP_KNL_LST','_EXP_KNL_DELTA_LST'};
rslts=cellfun(@(x)['MCMC_' IPD x '_SPRT_PARAS'],Mdls,'UniformOutput',false);
u=[{3,[1,3]},repmat({1:3,1:4},1,8)];
typ=[{'Const','Const'},repmat({'Cauchy','Cauchy','Exp','Exp'},1,2),{'Exp','Exp','Exp','Exp'},{'Cauchy','Cauchy','Exp','Exp'}];
beta0s=[0,beta0*ones(1,numel(rslts)-1)];
alpha0s=[0,0,alpha0*ones(1,numel(rslts)-2)];
delta0s=[0,0,repmat([0,delta0],1,8)];
h0s=[h0*ones(1,6),zeros(1,4),0.1,0.1,0.5,0.5,h0*ones(1,4)];
inclLST=[false(1,14),true(1,4)];

%% Plot pairwise distance distributions of all individuals and cases
PlotIndvdlDistDistns(db)

%% Run MCMC
RunVLStmMCMCSprtParas(db,r,a,b,u,beta0s,alpha0s,epsilon0,delta0s,h0s,typ,inclLST,niters,plotOutpt,rslts)

%% Process output
burnin=round(niters/10);
str='DIC'; % string to prepend to rslts when saving DIC
sprtParas='SprtParas';

%% Calculate DIC
RunCalcDICSprtParas(rslts(1:14),burnin,str)

%% Output results
RunCalcParEstsAndDICdiffs(rslts,burnin,IPD,str,false,false,sprtParas)

%% Plot deviance distributions
RunPlotDevDistn(rslts(1:6),burnin,IPD,sprtParas)