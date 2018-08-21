function RunVLStmMCMCSprtParas(db,r,a,b,u,beta0s,alpha0s,epsilon0,delta0s,h0s,typ,inclLST,niters,plotOutpt,rslts)
close all
%% RUN MCMC FOR DIFFERENT MODELS IN PARALLEL
nMdls=numel(rslts);

parfor i=1:nMdls
    VLStmMCMCSprtParas(db,r,a,b,u{i},beta0s(i),alpha0s(i),epsilon0,delta0s(i),h0s(i),typ{i},inclLST(i),niters,plotOutpt,rslts{i}) 
end