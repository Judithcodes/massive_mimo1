clear;clc;
M=20;
K=10;
p_db=0;
m=4;

N=10000;
pu=10^(p_db/10);
hMod = comm.RectangularQAMModulator('ModulationOrder',m,'BitInput',true,'NormalizationMethod','Average power');
hDemod = comm.RectangularQAMDemodulator('ModulationOrder',m,'BitOutput',true,'NormalizationMethod','Average power');
S=zeros(1,N*K);
R=zeros(1,N*K);
g1_rnd=zeros(1,N*K); %
g2_rnd=zeros(1,N*K); % complete formula
for i=0:N-1
    u=randi([0 1],K*log2(m),1);
    s=step(hMod,u);
    H=1/sqrt(2)*(randn(M,K)+1j*randn(M,K));
    n=1/sqrt(2)*(randn(M,1)+1j*randn(M,1));
    y=sqrt(pu)*H*s;
    A=H;%*inv(H'*H); %MRC
    r=(A'*y);
    xs=(real(s).*real(r)+imag(s).*imag(r))./(real(s).*real(s)+imag(s).*imag(s));
    xn=r-xs.*s;
    g1_rnd(i*K+1:i*K+K)=xs/K;
    AH=A'*H;
    num=transpose(pu*diag(AH).^2);
    den=pu*sum(abs(AH-diag(diag(AH))).^2)+sum(abs(H).^2);
    g2_rnd(i*K+1:i*K+K)=num./den;
        
%     S(i*K+1:i*K+K)=s;
%     R(i*K+1:i*K+K)=r;
    
end
num=pu*gamrnd(M,1,N,1);
den=1+pu*gamrnd(K-1,1,N,1);
g3_rnd=num./den;
[f3,g3]=hist(g3_rnd,100);
f2=(f3/length(g3_rnd)/(g3(2)-g3(1)));
bar(g3,f2);
hold on
[f2,g2]=hist(g2_rnd,100);
f2=f2/N/K/(g2(2)-g2(1));
plot(g2,f2);
[f1,g1]=hist(g1_rnd,100);
f1=f1/N/K/(g1(2)-g1(1));
plot(g1,f1,'linewidth',2);

MEAN=M*(exp(1/pu)*double(vpa(expint(K-1,sym(1/pu)),40)));
VAR=(M*(M+1))*((exp(1/pu)*(K-2+1/pu)*double(vpa(expint(K-2,sym(1/pu)),40))-1)/(K-2))-(MEAN)^2;
g=0:.1:10;
ft=-(g-1)/(2*pu)+log(g)*(M-1)-log(gamma(M))-log(pu)*((M+K-2)/2)-log(g+1)*((M+K)/2)+log(whittakerW((M-K+2)/2,(1-M-K)/2,(g+1)/pu));
plot(g,exp(ft),'k')
plot(g,gampdf(g,MEAN^2/VAR,VAR/MEAN),'b','linewidth',2);