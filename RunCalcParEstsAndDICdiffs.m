function RunCalcParEstsAndDICdiffs(rslts,burnin,IPD,str,doPlots,savePlots,sprtParas)
% set(0,'DefaultFigureVisible','off')
nMdls=numel(rslts);

%% Calculate parameter estimates and CIs for different models
thin=1; % no thinning currently to match DIC calcs
np=4;
mode_p=NaN(nMdls,np);
HPDI=NaN(np,2,nMdls);
mode_p1=NaN(nMdls,1);
HPDI1=NaN(nMdls,2);
pcorr=cell(nMdls,1);
mode_sptl=NaN(nMdls,1);
HPDI_sptl=NaN(nMdls,2);
mode_bckgrnd=NaN(nMdls,1);
HPDI_bckgrnd=NaN(nMdls,2);
mode_d_half=NaN(nMdls,1);
HPDI_d_half=NaN(nMdls,2);
mode_WHHRI=NaN(nMdls,1);
HPDI_WHHRI=NaN(nMdls,2);
mean_IP=NaN(nMdls,1);
HPDI_IP=NaN(nMdls,2);
MdlParEsts=NaN(nMdls,3*(np+1));
for i=1:nMdls
    fprintf([rslts{i} '\n'])
    [mode_p(i,:),HPDI(:,:,i),mode_p1(i),HPDI1(i,:),pcorr{i},mode_sptl(i),HPDI_sptl(i,:),mode_bckgrnd(i),HPDI_bckgrnd(i,:),mode_d_half(i),HPDI_d_half(i,:),mode_WHHRI(i),HPDI_WHHRI(i,:),mean_IP(i),HPDI_IP(i,:)]=ProcessOutput(rslts{i},burnin,thin,doPlots,savePlots,sprtParas);
    tmp=[];
    for j=1:np
        tmp=[tmp,mode_p(i,j),HPDI(j,:,i)];
    end
    MdlParEsts(i,:)=[tmp,mode_p1(i),HPDI1(i,:)];
end

save(['MdlParEstsFinal' IPD sprtParas])

%% Calculate DIC differences from best-fitting model
idx=1:6;
DICrslts=cellfun(@(x)[str '_' x],rslts(idx),'UniformOutput',false);
DICs=[];
for i=1:numel(DICrslts)
[DICs,DICmin,DICdiffs,RMLs]=CalcDICdiffs(DICrslts{i},DICs);
end

%% Output parameter estimates and DICs to file
MdlParEstsAndDICs=[MdlParEsts([idx,end-3,end-1],:),[DICs;NaN(2,1)],[DICdiffs;NaN(2,1)]];
save(['MdlParEstsAndDICs' IPD sprtParas],'MdlParEstsAndDICs')
ord=[1:3,5,4,6];
PrintModesAndCIsToFile(MdlParEstsAndDICs([ord,end-1,end],[1:end-5,end-1]),['MdlParEstsAndDICs' IPD sprtParas '.txt'])

%% Output acceptance rates
acc_rate_p=NaN(nMdls,1);
acc_rate_E=NaN(nMdls,1);
acc_rate_I=NaN(nMdls,1);
acc_rate_ERmove=NaN(nMdls,1);
acc_rate_R=NaN(nMdls,1);
for i=1:nMdls
    [acc_rate_p(i),acc_rate_E(i),acc_rate_I(i),acc_rate_ERmove(i),acc_rate_R(i)]=GetAccRates(rslts{i}); 
end
acc_rate=[acc_rate_p,acc_rate_E,acc_rate_I,acc_rate_ERmove,acc_rate_R];
save(['AccRates' IPD sprtParas],'acc_rate')

%% Presymptomatic infectiousness sensitivity analysis
%% Calculate DIC differences
idx1=[5:6,9:14];
DICrslts1=cellfun(@(x)[str '_' x],rslts(idx1),'UniformOutput',false);
DICs1=[];
for i=1:numel(DICrslts1)
[DICs1,DICmin1,DICdiffs1,RMLs1]=CalcDICdiffs(DICrslts1{i},DICs1);
end

%% Output parameter estimates and DICs to file
h0s=NaN(numel(idx1),1);
for i=1:numel(idx1)
    x=load(rslts{idx1(i)});
    h0s(i)=x.h0;
end
%%
MdlParEstsAndDICs1=[MdlParEsts(idx1,:),DICs1,DICdiffs1];
save(['MdlParEstsAndDICsh0SnstvtyAnlyss' IPD sprtParas],'MdlParEstsAndDICs1')
ord1=[3,1,5:2:numel(idx1),4,2,6:2:numel(idx1)];
PrintModesAndCIsToFile1([h0s(ord1),MdlParEstsAndDICs1(ord1,[1:end-5,end-1])],['MdlParEstsAndDICsh0SnstvtyAnlyss' IPD sprtParas '.txt'])