function plotFigure()
% 收敛精度与迭代次数的关系
close all ;
load jingdu.mat
iter = find(kll>0) ;
figure1 = figure;
axes1 = axes('Parent',figure1,...    
    'XTickLabel',{'1','10','100'},...
    'XScale','log',...
    'XMinorTick','on',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');
X1 = [iter', iter', iter'] ;
YMatrix1 = [kll(iter)', jc(iter)', mcc2(iter)'] ;
loglog1 = semilogx(X1,YMatrix1,'Parent',axes1);
set(loglog1(1),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','+','LineStyle','--','DisplayName','Kuhn et al.');
set(loglog1(2),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','*','LineStyle','--','DisplayName','Chae J');
set(loglog1(3),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','o','DisplayName','Zheng et al.');
xlabel('Iteration Number','FontWeight','bold','FontName','Times New Roman');
ylabel('Standard Error(\epsilon)','FontWeight','bold','FontName','Times New Roman');
legend(axes1,'show');
% 

% 收敛精度与噪声的关系
% 高斯噪声
close all ;
load zaosheng.mat
figure1 = figure;
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'0.001','0.010','0.100'},...
    'XScale','log',...
    'XMinorTick','on',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Times New Roman');
ylim([0.03,0.072])
box(axes1,'on');
hold(axes1,'all');
X1 = [zaos', zaos', zaos'] ;
X1(end,:) = [] ;
YMatrix1 = [kll', jc',mcc2'] ;
YMatrix1(end,:) = [] ;
loglog1 = loglog(X1,YMatrix1,'Parent',axes1);
set(loglog1(1),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','+','LineStyle','--','DisplayName','Kuhn et al.');
set(loglog1(2),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','*','LineStyle','--','DisplayName','Chae J');
set(loglog1(3),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','o','DisplayName','Zheng et al.');
xlabel('Standard Noise','FontWeight','bold','FontName','Times New Roman');
ylabel('Standard Error(\epsilon)','FontWeight','bold','FontName','Times New Roman');
legend(axes1,'show');


% 收敛精度与噪声的关系
% 椒盐噪声
close all ;
load zaosheng_nogauss.mat
figure1 = figure;
axes1 = axes('Parent',figure1,'YTickLabel',{'0.01','0.10','1.00'},...
    'YScale','log',...
    'YMinorTick','on',...
    'XTickLabel',{'0.001','0.010','0.100'},...
    'XScale','log',...
    'XMinorTick','on',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'all');
X1 = [zaos', zaos', zaos'] ;
YMatrix1 = [kll', jc', mcc2'] ;
loglog1 = loglog(X1,YMatrix1,'Parent',axes1);
set(loglog1(1),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','+','LineStyle','--','DisplayName','Kuhn et al.');
set(loglog1(2),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','*','LineStyle','--','DisplayName','Chae J');
set(loglog1(3),'Color','black','LineWidth',1,'MarkerSize',5,'Marker','o','DisplayName','Zheng et al.');
xlabel('Standard Noise','FontWeight','bold','FontName','Times New Roman');
ylabel('Standard Error(\epsilon)','FontWeight','bold','FontName','Times New Roman');
legend(axes1,'show');
end