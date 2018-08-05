% clear
% clc
% load('SliptHalf_pca_Cognition.mat')
% headers=header_Cog;
% refer to 'SubAnalysis.m' and 'SplitHalf_CV.m'

PaperDim=get(gcf,'PaperSize');
set(gcf,'PaperPosition',[0 0 PaperDim])
% Arbitrary Text, on background invisible axis
ht=axes('position',[0 0	1 1],'Visible','off');
text(ht,.3,.98,'Personality','FontSize',20)
text(ht,.6,.98, '65 variables','FontSize',14)

text(ht,.1,.95,'Factor Definition','FontSize',14)
% text(ht,.1,.3,'Missingness Report','FontSize',14)

% Plotting axis
hf=axes('position',[.14 .71 .7 .2]);
dom_n=sum(null_eigen_spec); % sum of null eigenspetrum
dom=sum(eigen_spec); % sum of eigen-spectrum
pca_n=0;
pca=0;
p=0;
for i=1:size(eigen_spec,2)
    pca_n=null_eigen_spec(i)/dom_n*100;
    plot(i,pca_n,'go')
    hold on
    pca=eigen_spec(i)/dom*100;
    plot(i,pca,'bo')
    hold on
    p=p+pca;
    plot(i,p,'ro')
    hold on
end
xlabel('Number of Principal Components (PCs)');
xlim([1 size(dd,2)])
%xticks([1:size(dd,2)])
ylabel('% of variance explained by each PC', 'Color','b');%# Add a label to the left y axis
ylim([0,100])
legend('Null eigen-spectrum', 'Covariance eigen-spectrum','Cumulative variance explained','Location','best')
yyaxis right
ylim([0,100])
ylabel('% of cumulative variance explained', 'Color', 'r')

title('Scree plot & Null eigen-spectrum')

%------------------------------------
hf=axes('position',[.14 .44 .7 .2]);
dom_CV=sum(mean(eigen_spec_CV));
eigen_spec_mean=mean(eigen_spec_CV);
p=0;
for i=1:size(eigen_spec_mean,2)
    pca=eigen_spec_mean(i)/dom_CV*100;
    p=p+pca;
    plot(i,p,'ro')
    hold on
    plot(rbar*100,'bo')
end
xlim([1 size(eigen_spec_mean,2)])
ylim([10 100])
%xticks([1:size(dd,2)])
legend({'In sample', 'Out sample'}, 'Location','SouthEast')
xlabel('Number of Principal Components (PCs)')
ylabel(gca,'% of cumulative variance explained');%# Add a label to the left y axis
title('PCs in 5-fold cross validation')

%------------------------------------
text(ht,.1,.38,'PC summary table','FontSize',14)
% Table
uitable('Data', [1 eigen_spec_mean(1)/sum(eigen_spec_mean)*100 rbar(1)*100;2 eigen_spec_mean(2)/sum(eigen_spec_mean)*100 ...
        (rbar(2)-rbar(1))*100; 3 eigen_spec_mean(3)/sum(eigen_spec_mean)*100 (rbar(3)-rbar(2))*100],...
        'RowName',[],...
        'ColumnName', {'PC', 'Train %var', 'Test %var'},...
        'Units', 'Normalized',...
        'Position',[0.09 0.27 .39 .09]);
% this table needs to be filled by hand
uitable('Data', [50 10 11;70 23 28; 90 44 48],...
        'RowName',[],...
        'ColumnName', {'>%var', '#PCs Train', '#PCs Test'},...
        'Units', 'Normalized',...
        'Position',[0.48 0.27 .39 .09]);    
      
% [bla,ind]=sort(abs([header_Alc{:,2}]), 'descend');
% weight=header_Alc(ind,:);

% uitable('Data', [weight(1,1) weight(2,1) weight(3,1);weight(1,2), weight(2,2) weight(3,2);...
%                  weight(4,1) weight(5,1) weight(6,1);weight(4,2), weight(5,2) weight(6,2);],...
%         'RowName',{'Ordered var name', 'Weights','Ordered var name', 'Weights'},...
%         'ColumnName', [],...
%         'Units', 'Normalized',...
%         'ColumnWidth', {107  107  107  107},...
%         'Position',[0.11 0.38 .8 .095]);   
       
         
print -dpng Report-Per.png
