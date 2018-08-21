function RunCalcDICSprtParas(rslts,burnin,str)
%% Calculate DICs for models in parallel
nMdls=numel(rslts);
parfor i=1:nMdls
    CalcDICSprtParas(rslts{i},burnin,'mode',str);
end