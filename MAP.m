function M=MAP(p)
if all(p==p(1))
    M=p(1);
else
    [N,E]=histcounts(p,50);
    [~,i]=max(N);
    M=(E(i)+E(i+1))/2;
end