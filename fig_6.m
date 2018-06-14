clear;clc;
load('ber_resultsK50M400_200_100.mat')
M=[400 200 100];
K=50;
marker='^so';
for m=1:3
    semilogy(p_db(mod(p_db,2)==0),ber_sim(mod(p_db,2)==0,m),[marker(m) 'k'],'markersize',10,'markerfacecolor','k');
    hold on;
    semilogy(p_db,ber_the(:,m),'-b','linewidth',2');
    semilogy(p_db,ber_app(:,m),'--r','linewidth',2');
end
axis([-20 10 1e-3 1])
set(gca,'fontsize',18);
xlabel('$p_u$(dB)','fontsize',18,'interpreter','latex')
ylabel('BER','fontsize',18,'interpreter','latex')
legend('Simulation','Theory','ApproximationM=100');
dim = [.1 .8 .1 .1];
str = {[native2unicode(char(hex2dec('25cf')),'UTF-16') '  M=' num2str(100) ', K=50'],...
       [native2unicode(char(hex2dec('25a0')),'UTF-16') '  M=' num2str(200) ', K=50'],...
       [native2unicode(char(hex2dec('25b2')),'UTF-16') '  M=' num2str(400) ', K=50']};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',16,'BackgroundColor','w');
grid;
matlab2tikz('D:\OneDrive\Research\Massive Mimo\_Papers\1\New folder\IEEEtran\fig6.tex','width','\figW','height','\figH');

% clear
% clc
% antenna=[100 200 400];
% k=20;
% m=4;
% p_db=-22:2:-6;
% %p_db=10.^(p_b./10)*log2(m);
% n=10^5;
% AP_MRC_CLOSE=zeros(length(antenna),length(p_db));
% AP_MRC=zeros(length(antenna),length(p_db));
% BER_MRC=zeros(length(antenna),length(p_db));
% for j=1:length(antenna)
%     M=antenna(j);
%     K=k;
%     for i=1:(length(p_db))
%         p=10^(p_db(i)/10);
%         MEAN=M*(exp(1/p)*double(vpa(expint(K-1,sym(1/p)),40)));
%         m_2=(exp(1/p)*(K-2+1/p)*double(vpa(expint(K-2,sym(1/p)),40))-1)/(K-2);
%         v_m=(M*(M+1))*(m_2)-(MEAN)^2;
%         beta=v_m/MEAN;
%         alpha=MEAN/beta;
%         
%         SINR=gamrnd(M,1,n,1)./(gamrnd(K-1,1,n,1)+1/p);
%         BER_MRC(i,j)=(4*(sqrt(m)-1))/(log2(m)*(sqrt(m)))*mean(qfunc(sqrt((3*SINR)/(m-1))));
%         SINR=gamrnd(alpha,beta,n,1);
%         AP_MRC(i,j)=(4*(sqrt(m)-1))/(log2(m)*(sqrt(m)))*mean(qfunc(sqrt((3*SINR)/(m-1))));
%         AP_MRC_CLOSE(i,j)=.2/(1.5*beta/(m-1)+1)^(alpha);
%     end
% end
% 
%  %% plot 
% figure
% hold on
% f1=semilogy(p_db,BER_MRC,'r-s','LineWidth',2,'MarkerSize',10);
% f2=semilogy(p_db,AP_MRC,'b:*','LineWidth',2,'MarkerSize',10);
% f3=semilogy(p_db, AP_MRC_CLOSE,'b:>','LineWidth',2,'MarkerSize',10);
% ax = gca; % current axes
% ax.FontSize = 18;
% f1_group= hggroup;
% %f2_group= hggroup;
% f3_group= hggroup;
% set(f1,'Parent',f1_group)
% %set(f2,'Parent',f2_group)
% set(f3,'Parent',f3_group)
% % Include these hggroups in the legend:
% set(get(get(f1_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% %set(get(get(f2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% set(get(get(f3_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
% grid on
% box on
% xlabel('p_u (dB)','fontsize',18)
% ylabel(' BER ','fontsize',18)
% axis([min(p_db) max(p_db) 10^-5 1])
% legend('Simulation','Theory')
% dim = [0.2 0.2 0.25 0.25];
% str = {'MRC Receiver','K=10','4-QAM'};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',18,'LineWidth',2,'BackgroundColor','w');
% 
% % X = [.5 0.5];Y = [.7 .2];
% % txtar = annotation('textarrow',X,Y,'String',' K=5 ','FontSize',18,'LineWidth',2,'TextBackgroundColor','w','TextEdgeColor','K');
% % dim = [.5 0.15 .02 .05];
% % annotation('ellipse',dim,'LineWidth',2);
% 
% X = [0.6 0.6];Y = [.7 .31];
% txtar = annotation('textarrow',X,Y,'String',' M=400 ','FontSize',18,'LineWidth',2,'TextBackgroundColor','w','TextEdgeColor','K');
% dim = [.59 0.25 .02 .065];
% annotation('ellipse',dim,'LineWidth',1);
% 
% X = [0.7 0.7];Y = [.7 .48];
% annotation('textarrow',X,Y,'String',' M=200 ','FontSize',18,'LineWidth',2,'TextBackgroundColor','w','TextEdgeColor','K');
% dim = [.69 0.42 .02 .065];
% annotation('ellipse',dim,'LineWidth',1);
% 
% X = [0.8 0.8];Y = [.7 .62]; annotation('textarrow',X,Y,'String',' M=100 ','FontSize',18,'LineWidth',2,'TextBackgroundColor','w','TextEdgeColor','K');
% dim = [.79 0.56 .02 .065];
% annotation('ellipse',dim,'LineWidth',1);