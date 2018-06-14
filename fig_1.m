clear
clc
M=10;
K=5;
n=10^6;
p_db=[-10,0,10,20];


%% SIMULATION
for i=1:length(p_db)
    x=linspace(0,3*i,300);
   p=10^(p_db(i)/10); 
 num=p*gamrnd(M,1,n,1);
 den=1+p*gamrnd(K-1,1,n,1);
 Dd=num./den;
[count1,bins1]=hist(Dd,50*i^2);
f2=(count1/length(Dd)/(bins1(2)-bins1(1)));


%% TEORY
 PDF=zeros(1,length(x));
  PDF2=zeros(1,length(x));
for j=1:length(x)
    g=x(j);
    u=1/p:.01:100;
    y=u.^(M).*(u-1/p).^(K-2).*exp(-(g+1)*u);
    out1=sum(y)*(u(2)-u(1));
    PDF(j)=exp(1/p)*g^(M-1)/(gamma(M)*gamma(K-1))*out1;
    PDF2(j)=gamma(M+K-1)/gamma(M)/gamma(K-1)*g^(M-1)/(1+g)^(M+K-1);
end 

%%
subplot(2,length(p_db)/2,i)
hold on
  F2=bar(bins1,f2,'FaceColor','w','EdgeColor','k','LineWidth',2);
  f2=semilogx(x,PDF,'b','LineWidth',2);
  semilogx(x,PDF2,'--r','LineWidth',2);
  %f3=semilogx(x,aproximat,'k:','LineWidth',3,'MarkerSize',10')
  ax = gca; % current axes
  grid on
  box on
  xlabel('$\gamma$','interpreter','latex')
  ylabel('$f_{\gamma_k}(\gamma)$','interpreter','latex')
   %legend('Simulation','Theory ','Approximation ')
   legend('Simulation','Theory','Approx. in [16]')
axis([0 max(x) 0 max(PDF)+.1])
end
dim = [0.345 0.76 0.07 0.12];
str = {'p_u=-10 dB ','M=10, K=5' };
annotation('textbox',dim,'String',str,'FitBoxToText','off','LineWidth',2,'BackgroundColor','w');

dim = [0.786 0.76 0.07 0.12];
str = {'p_u=0 dB ','M=10, K=5' };
annotation('textbox',dim,'String',str,'FitBoxToText','off','LineWidth',2,'BackgroundColor','w');

dim = [0.345 0.29 .07 .12];
str = {'p_u=10 dB ','M=10, K=5' };
annotation('textbox',dim,'String',str,'FitBoxToText','off','LineWidth',2,'BackgroundColor','w');

dim = [0.786 0.29 0.07 0.12];
str = {'p_u=20 dB ','M=10, K=5' };
% annotation('textbox',dim,'String',str,'FitBoxToText','off','LineWidth',2,'BackgroundColor','w');
% matlab2tikz('D:\OneDrive\Research\Massive Mimo\_Papers\1\New folder\IEEEtran\fig2.tex','width','\figW','height','\figH');
close all
