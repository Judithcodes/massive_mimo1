clear
clc
M=100;
K=40;
n=10^6;
p_db=[-10,0,10,20];
%% SIMULATION
for i=1:length(p_db)
    x=linspace(0,6.5,40);
   p=10^(p_db(i)/10); 
 num=p*gamrnd(M,1,n,1);
 den=1+p*gamrnd(K-1,1,n,1);
 Dd=num./den;
[count1,bins1]=hist(Dd,25+5*i);
f2=(count1/length(Dd)/(bins1(2)-bins1(1)));


%% TEORY
PDF=zeros(1,length(x));
for j=1:length(x)
    g=x(j);
    u=1/p:.01:100;
    y=u.^(M).*(u-1/p).^(K-2).*exp(-(g+1)*u);
    out1=sum(y)*(u(2)-u(1));
    PDF(j)=exp(1/p)*g^(M-1)/(gamma(M)*gamma(K-1))*out1;
end 
%% GAMMA
m_1=(exp(1/p)*mfun('Ei',K-1,1/p));
m=M*m_1;
m_2=(exp(1/p)*(K-2+1/p)*mfun('Ei',K-2,1/p)-1)/(K-2);
v_m=(M*(M+1))*(m_2)-((M)*(m_1))^2;
Beta=v_m/m;
alpha=m/Beta;
aproximat=gampdf(x,alpha,Beta);

%% UPPER BOUNDS
UPPER_BOUN=1/beta(M,K-1)*x.^(M-1)./(1+x).^(M+K-1);

%%
subplot(2,length(p_db)/2,i)
hold on
  bar(bins1,f2,'FaceColor','w','EdgeColor','k','LineWidth',2)
  semilogx(x,PDF,'r-*','LineWidth',2,'MarkerSize',10')
  semilogx(x,aproximat,'b-S','LineWidth',2,'MarkerSize',10')
  ax = gca; % current axes
  ax.FontSize = 18;
  grid on
  box on
  xlabel('\gamma_k','fontsize',18)
  ylabel('Probability Density','fontsize',18)
  legend('Simulation','Theory ','Approximation ')
axis([0 max(x) 0 max(PDF)+.1])
end
dim = [0.32 0.30 0.5 0.45];
str = {'p_u=-10 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');

dim = [0.77 0.66 0.1 0.1];
str = {'p_u=0 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');

dim = [0.32 0.28 .01 .01];
str = {'p_u=10 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');

dim = [0.77 0.10 0.20 0.19];
str = {'p_u=20 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');


