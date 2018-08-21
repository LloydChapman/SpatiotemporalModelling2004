function RunPlotDevDistn(rslts,burnin,IPD,str)
figure;
nMdls=numel(rslts);

for i=1:nMdls
    hL(i)=PlotDevDistn(rslts{i},burnin); hold on
end
legend(hL,cellfun(@(x)x(11:end),rslts,'UniformOutput',false),'Interpreter','none')
saveas(gcf,['PstrDevDistns' IPD str])
saveaspdf(gcf,['PstrDevDistns' IPD str])