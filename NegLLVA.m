function NLL=NegLLVA(pars,data,cens,freq)
%NEGLL_VA Negative log-likelihood function for variable asymptote catalytic model

% Reshape data input from vector to matrix
data=[data(1:4:end),data(2:4:end),data(3:4:end),data(4:4:end)]';

lambda=pars(1); c=pars(2);
n=data(3,:); 
k=data(4,:); 
a=(data(1,:)+data(2,:))/2;
p=PrevVA(a,lambda,c);
LL_a=k.*log(p)+(n-k).*log(1-p);
NLL=-sum(LL_a);