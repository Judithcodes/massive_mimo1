%test_expint
clear;
clc;
p=1;
K=50;

i1=eval(vpa(exp(1/p)*expint(K-1,sym(1/p)),25));

n=K-2;
b=1/p;
i2=(-1)^n*b^n*exp(b)*expint(b);
for k=1:n
    i2=i2+factorial(k-1)*(-b)^(n-k);
end;
i2=i2/factorial(K-2);
fprintf('%f - %f  = %e\n',i1,i2,i1-i2);