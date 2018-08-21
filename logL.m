function [LL,LL1,LL2,LL3]=logL(S,lambda,tEm,IP,r,p)
probS=S.*lambda;
LL1=-sum(probS(:)); % probabilities of staying susceptible
LL2=sum(log(1-exp(-lambda(tEm)))); % probabilities of being infected in given months
LL3=sum(log(nbinpdf(IP-1,r,p))); % probabilities of incubation period durations
LL=LL1+LL2+LL3;
