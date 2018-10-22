function [M,HPDI]=PlotPstrDistn(p,pname,nbins,prior,varargin)
h=histogram(p,nbins,'Normalization','pdf'); hold on
[y,i]=max(h.Values);
if all(p==p(1))
    M=p(1);
    HPDI=[p(1) p(1)];
else
    M=(h.BinEdges(i)+h.BinEdges(i+1))/2;
    probs=h.Values.*h.BinWidth;
    sort_probs=sort(probs,'descend');
    HPDIht_idx=find(cumsum(sort_probs)>=0.95,1);
    HPDIht=sort_probs(HPDIht_idx);
    indcs=find(probs>=HPDIht);
    HPDI=[(h.BinEdges(indcs(1))+h.BinEdges(indcs(1)+1))/2 (h.BinEdges(indcs(end))+h.BinEdges(indcs(end)+1))/2];
end
% Plot HPDI
plot([HPDI(1) HPDI(1)],[0 1.05*y],'r',[HPDI(2) HPDI(2)],[0 1.05*y],'r','LineWidth',1.5)
% indcs=min(indcs):max(indcs); % make sure range for HPDI is continuous
% bar((hp.BinEdges(indcs)+hp.BinEdges(indcs+1))/2,hp.Values(indcs),1,'FaceColor',[0.85 0.325 0.098])
% Plot mode
plot([M M],[0 1.05*y],'m','LineWidth',1.5)
% Plot prior distn
if nargin==5
    a=varargin{1};
    plot(h.BinEdges,pdf(prior,h.BinEdges,a),'g','LineWidth',1.5)
elseif nargin==6
    a=varargin{1};
    b=varargin{2};
    plot(h.BinEdges,pdf(prior,h.BinEdges,a,b),'g','LineWidth',1.5)
elseif nargin==7
    a=varargin{1};
    b=varargin{2};
    c=varargin{3};
    plot(h.BinEdges,pdf(prior,h.BinEdges,a,b,c),'g','LineWidth',1.5)
end
axis([min(h.BinEdges) max(h.BinEdges) 0 1.05*y])
% set(gca,'Fontsize',16)
if ismember(pname,{'beta','alpha','epsilon','delta'})
    xlabel(['\' pname],'FontSize',13)
%     xlabel(['\' pname],'FontSize',20)
else
    xlabel(['$$' pname '$$'],'FontSize',13,'Interpreter','latex')
%     xlabel(pname,'FontSize',16);
end
ylabel('Density')
% ylabel('Probability density','FontSize',14)
hold off
% saveas(gcf,['PosteriorDistn' pname{i} '.eps'],'epsc')