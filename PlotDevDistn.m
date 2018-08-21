function h=PlotDevDistn(str,burnin1)
% Load MCMC output
load(str)

% Redefine burnin if necessary
if ~exist('z','var')
    z=burnin1+1:niters;
end

% DIC=dic(-2*LL(end-numel(LLM)+1:end),-2*LLM);
% histogram(DIC,'Normalization','pdf')

% Plot deviance distribution
h=histogram(-2*LL(z),50,'Normalization','pdf');
xlabel('Posterior deviance'); ylabel('Density')
