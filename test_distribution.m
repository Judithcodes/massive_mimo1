clear
clc
M=20;
K=5;
n=10^6;
p_db=-40;
p=10^(p_db/10);

num=p*gamrnd(M,1,n,1);
den=1+p*gamrnd(K-1,1,n,1);
Dd=num./den;
[count1,bins1]=hist(Dd,50);
f2=(count1/length(Dd)/(bins1(2)-bins1(1)));
x=linspace(0,max(bins1),400);

PDF=(-1)^M*M*exp(-x/p).*x.^(M-1)./(x+1).^(M+K-1).*laguerreL(M,1-K-M,(1+x)/p);
PDF2=exp(-(x-1)/2/p).*x.^(M-1)/gamma(M)/p^((M+K-2)/2)./(x+1).^((M+K)/2).*whittakerW((M-K+2)/2,(1-M-K)/2,(x+1)/p);
aproximat2=M*exp(-x/p).*x.^(M-1)./(x+1).^(M+K-1)/gamma(M+1).*(x+1).^(M+K-1)/p^M;
aproximat3=M*exp(-x/p).*x.^(M-1)./(x+1).^(M+K-1)/gamma(M+1).*(x+1).^(M)/p^M;
m_1=(exp(1/p)*double(vpa(expint(sym(K-1),1/p))));
m=M*m_1;
mmm=M/(K-2+1/p);
m_2=(exp(1/p)*(K-2+1/p)*double(vpa(expint(sym(K-2),1/p)))-1)/(K-2);
v_m=(M*(M+1))*(m_2)-((M)*(m_1))^2;
vvv=M*(M+1)*(1/(K-3+1/p)-1/(K-2+1/p))-mmm^2;
Beta=v_m/m;
Alpha=m/Beta;
bbb=vvv/mmm;
aaa=mmm/bbb;
aproximat1=gampdf(x,aaa,bbb);
% aproximat2=exp(-x/p).*x.^(M-1)/gamma(M)/p^M;%gampdf(x,M,p);

hold on
bar(bins1,f2,'FaceColor','w','EdgeColor','k','LineWidth',1)
semilogx(x,PDF,'--r','LineWidth',3)
semilogx(x,PDF2,'g','LineWidth',5)
semilogx(x,aproximat1,'b','LineWidth',3)
semilogx(x,aproximat2,'--k','LineWidth',3)
semilogx(x,aproximat2,'y','LineWidth',2)

grid on
box on
xlabel('$\gamma$','interpreter','latex')
ylabel('$f_{\gamma_k}(\gamma)$','interpreter','latex')
legend('Simulation','Theory ','Approximation ')
% figure
% plot(x,whittakerW((M-K+2)/2,(1-M-K)/2,(x+1)/p));
% hold on
% plot(x,exp(-(x+1)/2/p).*(1+x).^(-(M+K)/2)/p^((M-K+2)/2),'k');
% plot(x,exp(-(x+1)/2/p).*((x+1)/p).^((M-K+2)/2),'r');
figure
plot(x,laguerreL(M,1-K-M,(1+x)/p));
hold on
plot(x,(-1)^M/gamma(M+1).*(x+1).^(M+K-1)/p^M,'k');
plot(x,(-1)^M/gamma(M+1).*(x+1).^(M)/p^M,'r');


