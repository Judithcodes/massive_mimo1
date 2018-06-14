clear;clc;
% % a=-5;
% % n=50;
g=1.5;
p=.01;
z=(1+g)/p;
M=200;
K=10;


a=1-K-M;
n=M;
w1=(-1)^n*factorial(n)*exp(-1/2*z)*z^(a/2+.5)*laguerreL(n,a,z);
w2=whittakerW(a/2+.5+n,a/2,z);




% G=0:.01:10;
% pdf3=M*exp(-G/p).*G.^(M-1)./(G+1).^(M+K-1)/factorial(M).*((1+G)/p).^M;
% pdf1=zeros(size(G));
% pdf3=zeros(size(G));
% for i=1:length(G)
%     g=G(i);
    u=1/p:.01:1000;
    y=log(u)*(M)+log(u-1/p)*(K-2)+(-(g+1)*u)-sum(log(1:M-1))-sum(log(1:K-2));
    out1=sum(exp(y))*(u(2)-u(1));
    PDF1=exp(1/p)*g^(M-1)*out1;
    PDF2=exp(-(g-1)/2/p)*g^(M-1)/gamma(M)/p^((M+K-2)/2)/(g+1)^((M+K)/2)*w1;
    PDF3=(-1)^M*M*exp(-g/p)*g^(M-1)/(g+1)^(M+K-1)*laguerreL(n,a,z);
    
%     PDF4=    M*exp(-g/p)*g^(M-1)/(g+1)^(M+K-1)/factorial(M)*((1+g)/p)^M;
%     PDF5= exp((-g/p)+log(g)*(M-1)-log(g+1)*(K-1)-sum(log(1:M-1))-log(p)*M);
%     pdf1(i)=PDF1;
% end;