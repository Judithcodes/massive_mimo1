clear
clc
M=100;
K=40;
n=10^6;
p_db=[-10,0,10,20];
%% SIMULATION
for i=1:length(p_db)
    x=linspace(0,4.5,30);
    p=10^(p_db(i)/10);
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
    
    %%
    subplot(2,length(p_db)/2,i)
    hold on
    semilogx(x,p_outage,'sb','LineWidth',2,'MarkerSize',6)
    semilogx(x,AP_P_outage,'-b','LineWidth',2,'MarkerSize',6)
    ax = gca; % current axes
    ax.FontSize = 18;
    grid on
    box on
    xlabel('$\gamma_{th}$','fontsize',18,'interpreter','latex')
    ylabel('Outage Probability','fontsize',18,'interpreter','latex')
    legend('Simulation','Approximation','Location','southeast')
    axis([0 max(x) 0 1])
end
dim = [0.32 0.80 0.01 0.01];
str = {'p_u=-10 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');

dim = [0.77 0.80 0.01 0.01];
str = {'p_u=0 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');

dim = [0.32 0.5 .01 .01];
str = {'p_u=10 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');

dim = [0.77 0.5 0.01 0.01];
str = {'p_u=20 dB ',' M=100, K=40' };
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');
matlab2tikz('D:\OneDrive\Research\Massive Mimo\_Papers\1\New folder\IEEEtran\fig5.tex','width','\figW','height','\figH');
% close all

