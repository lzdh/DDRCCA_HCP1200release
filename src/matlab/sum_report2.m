% run 'SDRvariables.m'

load('workspace/SDRvariables.mat')

%sub_area={1Demo, 2PhyHeal, 3FemHeal, 4FamHist, 5Psychiatry, 6Sensory, 7Drug, 8Alcohol, 9Tobacco, 10Alertness, 11Cognition, 
%          12Emotion, 13Motor, 14Personality};
Error1 = Error1{14};
Error2 = Error2{14};
%sub_header = erase(sub_header,"FamHist_");


PaperDim=get(gcf,'PaperSize');
set(gcf,'PaperPosition',[0 0 PaperDim])
% Arbitrary Text, on background invisible axis
ht=axes('position',[0 0	1 1],'Visible','off');
text(ht,.4,.99,'Personality','FontSize',20)

text(ht,.1,.97,'Variable weights','FontSize',14)

hf=axes('position',[.06 .74 .89 .2]);
plot(weight(:,1:3),'linewidth',2);
hold on
stem(weight(:,1:3),':', 'linewidth',1)
set(gca,'XTickLabel',sub_header,'XTickLabelRotation',45,'XGrid','on',...
             'XTick',1:size(sub_header,1),'TickLabelInterpreter','none','FontName', 'Arial Narrow','Fontsize',9);
%ax = gca;
%ax.TickLabelInterpreter = 'latex';
%ax.XTickLabel{1} = ['\color{red}' ax.XTickLabel{1}];
xlim([0 size(dd,2)+1])
legend(strcat({'PC '},num2str([1:3]')),'Location','best','Fontsize',7) %[0.45 0.15 0.05 0.03]
ylabel('PC weight')
title('First 3 Pricipal Loadings')  



hf=axes('position',[.06 .44 .89 .19]);
plot(weight_obl,'linewidth',2);
hold on
stem(weight_obl(:,1:3),':', 'linewidth',1)
set(gca,'XTickLabel',sub_header(:,1),'XTickLabelRotation',45,'XGrid','on',...
             'XTick',1:size(sub_header,1),'TickLabelInterpreter','none','FontName', 'Arial Narrow','Fontsize',9);
legend(strcat({'RC '},num2str([1:3]')),'Location','south','Fontsize',7)
ylabel('RC weight')
xlim([0 size(dd,2)+1])
title('First 3 Rotated Component (RC) Loadings','Interpreter','none')


text(ht,.1,.32,'Optimal number of PCs','FontSize',14)

hf=axes('position',[.14 .10 .74 .19]);
plot(Error1(1:15), 'k.--')
hold on
plot(Error2(1:15), 'r.-')
legend({ 'Simple method', 'Two-way CV method'}, 'Location', 'best')
legend boxoff
set(gca, 'XTick', 1:15)
%set(gca, 'YTick', [])
xlim([1 15])
xlabel('Number of PCs')
ylabel('Cross-validation error')
title('Error analysis using different number of PCs') 


print -dpng Report-Per2.png  