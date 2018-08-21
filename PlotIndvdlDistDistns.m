function PlotIndvdlDistDistns(str)
load(str)
dists=CalcDists(data);
KAdists=dists(data.KA==1);
npara=numel(unique(data.Para));
figure;
PlotDistn(dists,KAdists,'All paras'); hold on

for i=1:npara
dists1=dists(data.Para==i,data.Para==i);
KAdists1=dists((data.Para==i&data.KA==1),(data.Para==i&data.KA==1));
figure;
PlotDistn(dists1,KAdists1,['Para ' num2str(i)])
end

function PlotDistn(x,y,str)
x(1:size(x,1)+1:end)=NaN;
y(1:size(x,1)+1:end)=NaN;
h=histogram(x,40,'Normalization','probability'); hold on
histogram(y,h.BinEdges,'Normalization','probability')
set(gca,'FontSize',20)
xlabel('Distance (m)'); ylabel('Density'); title(str,'Interpreter','none')
legend('All individuals','VL cases')
hold off
saveas(gcf,['INDVDL_dists_distn_' str])
saveaspdf(gcf,['INDVDL_dists_distn_' str])