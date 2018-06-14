clear
clc
MM=[50 50 200 200];
KK=[10 40  40 180];
max_x=[18 3.5 14 3];
max_y=[.3 2 .5 4];
nbins=[40 25 20 12];
n=10^6;
p_db=-20;
%% SIMULATION
for i=1:length(MM)

 p=10^(p_db/10);
 M=MM(i);
 K=KK(i);
 num=p*gamrnd(M,1,n,1);
 den=1+p*gamrnd(K-1,1,n,1);
 Dd=num./den;
[count1,bins1]=hist(Dd,nbins(i));
f2=(count1/length(Dd)/(bins1(2)-bins1(1)));
x=linspace(0,max_x(i),400);

%% TEORY
PDF=zeros(1,length(x));
PDF2=zeros(1,length(x));
for j=1:length(x)
    g=x(j);
    u=1/p:.01:10000;
    y=log(u)*(M)+log(u-1/p)*(K-2)+(-(g+1)*u)-sum(log(1:M-1))-sum(log(1:K-2));
    out1=sum(exp(y))*(u(2)-u(1));
    PDF(j)=exp(1/p)*g^(M-1)*out1;
    PDF2(j)=(-1)^M*M*exp(-g/p)*g^(M-1)/(g+1)^(M+K-1)*laguerreL(M,1-K-M,(1+g)/p);
end 
%% GAMMA
m_1=(exp(1/p)*double(vpa(expint(sym(K-1),1/p))));
m=M*m_1;
mmm=M/(K-2+1/p);
m_2=(exp(1/p)*(K-2+1/p)*double(vpa(expint(sym(K-2),1/p)))-1)/(K-2);
% mmm2=1/(K-3+1/p)-1/(K-2+1/p);
v_m=(M*(M+1))*(m_2)-((M)*(m_1))^2;
vvv=M*(M+1)*(1/(K-3+1/p)-1/(K-2+1/p))-mmm^2;
Beta=v_m/m;
Alpha=m/Beta;
bbb=vvv/mmm;
aaa=mmm/bbb;
aproximat=gampdf(x,aaa,bbb);

%% UPPER BOUNDS
UPPER_BOUN=1/beta(M,K-1)*x.^(M-1)./(1+x).^(M+K-1);

%%
subplot(2,length(MM)/2,i)
hold on
  bar(bins1,f2,'FaceColor','w','EdgeColor','k','LineWidth',1)
  semilogx(x,PDF,'--r','LineWidth',3)
  semilogx(x,PDF2,'--k','LineWidth',1)
  semilogx(x,aproximat,'b','LineWidth',3)
  ax = gca; % current axes
%   ax.FontSize = 14;
  grid on
  box on
  xlabel('$\gamma$','interpreter','latex')
  ylabel('$f_{\gamma_k}(\gamma)$','interpreter','latex')
  legend('Simulation','Theory ','Approximation ')
% axis([0 max_x(i) 0 max_y(i)])
end
dim = [0.34 0.76 0.07 0.12];
str = {'p_u=0 dB','M=50, K=10' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',2,'BackgroundColor','w');

dim = [0.78 0.76 0.07 0.12];
str = {'p_u=0 dB','M=50, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',2,'BackgroundColor','w');

dim = [0.33 0.28 .07 .12];
str = {'p_u=0 dB','M=200, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',2,'BackgroundColor','w');

dim = [0.76 0.28 0.07 0.12];
str = {'p_u=0 dB','M=200, K=180' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineWidth',2,'BackgroundColor','w');
% matlab2tikz('D:\OneDrive\Research\Massive Mimo\_Papers\1\New folder\IEEEtran\fig3.tex','width','\figW','height','\figH');
% close all


