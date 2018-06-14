clear
clc
M=100;
K=40;
n=10^6;
p_db=10;
%% SIMULATION
x=linspace(0,4.5,30);
p=10^(p_db/10);
AP_P_outage=zeros(1,length(x));
p_outage=zeros(1,length(x));
for j=1:length(x)
    g_T=x(j);
    num=p*gamrnd(M,1,n,1);
    den=1+p*gamrnd(K-1,1,n,1);
    SINR=num./den;
    sinr=SINR(SINR<g_T);
    p_outage(j)=length(sinr)/n;
    
    m_1=(exp(1/p)*double(vpa(expint(sym(K-1),1/p),40)));
    m=M*m_1;
    m_2=(exp(1/p)*(K-2+1/p)*double(vpa(expint(sym(K-2),1/p),40))-1)/(K-2);
    v_m=(M*(M+1))*(m_2)-((M)*(m_1))^2;
    Beta=v_m/m;
    alpha=m/Beta;
    AP_P_outage(j)=1-igamma(alpha,g_T/Beta)/gamma(alpha);
end
hold on
semilogx(x,p_outage,'sb','LineWidth',2,'MarkerSize',6)
semilogx(x,AP_P_outage,'-b','LineWidth',2,'MarkerSize',6)
ax = gca; % current axes
ax.FontSize = 18;
grid on
box on
xlabel('$\gamma_{th}$','fontsize',18,'interpreter','latex')
ylabel('Outage Probability','fontsize',18,'interpreter','latex')
legend('Simulation','Approximation','Location','northwest')
axis([0 max(x) 0 1])
amatlab2tikz('D:\OneDrive\Research\Massive Mimo\_Papers\1\New folder\IEEEtran\fig5.tex','width','\figW','height','\figH');
% close all

