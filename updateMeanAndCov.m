function [mean_new,cov_new]=updateMeanAndCov(mean_old,cov_old,p,i)
mean_new=(i-1)/i*mean_old+p/i;
cov_new=(i-1)/i*(cov_old+(p-mean_old)'*(p-mean_old)/i);