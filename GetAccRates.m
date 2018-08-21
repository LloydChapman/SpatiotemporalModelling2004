function [acc_rate_p,acc_rate_E,acc_rate_I,acc_rate_ERmove,acc_rate_R]=GetAccRates(data)
S=load(data);
acc_rate_p=S.acc_rate_p;
acc_rate_E=S.acc_rate_E;
acc_rate_I=S.acc_rate_I;
acc_rate_ERmove=S.acc_rate_ERmove;
acc_rate_R=S.acc_rate_R;