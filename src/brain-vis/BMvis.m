I=2;

% ts_dir=sprintf('/Volumes/Jessie/HCP/HCP_PTN820/node_timeseries/3T_HCP820_MSMAll_d200_ts2');
% ts = nets_load(ts_dir,0.72,1,1);

% to get the interactive connectome viewer (note this is bundled with fslnets now)
NM='/Users/jessie_liu/Desktop/HCP1200/BrainVis';
SUMPICS='/Users/jessie_liu/Desktop/HCP1200/rawdata/groupICA/melodic_IC_sum.sum'; 
%SUMPICS='/Volumes/Jessie/HCP/HCP_PTN820/groupICA/groupICA_3T_HCP820_MSMAll_d200.ica/melodic_IC_sum.sum'; % thicker versions of the thumbnails

netmat1=load('../rawdata/netmats1.txt');
netmat2=load('../rawdata/netmats2.txt');


% netmats1=load('/Volumes/Jessie/HCP/HCP_PTN820/netmats/3T_HCP820_MSMAll_d200_ts2/netmats1.txt');
% netmats2=load('/Volumes/Jessie/HCP/HCP_PTN820/netmats/3T_HCP820_MSMAll_d200_ts2/netmats2.txt');


[Znet1,Mnet1]=nets_groupmean(netmat1,0,1); % full correlation
[Znet2,Mnet2]=nets_groupmean(netmat2,0,1); % partial correlation

Znet1(Znet1 > max(Znet1(~isinf(Znet1))))=max(Znet1(~isinf(Znet1)));
Znet1(Znet1 < min(Znet1(~isinf(Znet1))))=min(Znet1(~isinf(Znet1)));

Znet2(Znet2 > max(Znet2(~isinf(Znet2))))=max(Znet2(~isinf(Znet2)));
Znet2(Znet2 < min(Znet2(~isinf(Znet2))))=min(Znet2(~isinf(Znet2)));

%heatmap on correlations
heatmap(Znet1,'Colormap',jet)
heatmap(Znet2,'Colormap',jet)


%%% all-200-nodes hierarchical plot
[hierALL,linkagesALL]=nets_hierarchy(Znet1,Znet2,[1:200],SUMPICS,1.5);
clustersALL=cluster(linkagesALL,'maxclust',4);
set(gcf,'PaperPositionMode','auto','Position',[1 1 46*(200+1) 1574 ]);


ZnetMOD=reshape(CL_bm(:,1),200,200);
[hierALL,linkagesALL]=nets_hierarchy(Znet1, ZnetMOD*300, [1:200],SUMPICS,1.5);

%heatmap(ZnetMOD,'Colormap',jet)


%%% select top-30 CCA-edges (list of nodes goes into grotDD)
grot=ZnetMOD; 
grotTHRESH=prctile(abs(grot(:)),99.5); % 99.85% is about 30 edges 19800*(1-0.9985); 99.5% is about 100 edges
grot(abs(grot)<grotTHRESH)=0;  grot(isnan(grot))=0;
grotDD=find(sum(grot~=0)>0);  grot=grot(grotDD,grotDD); grotTHRESH
grot1=Znet1; grot1=grot1(grotDD,grotDD); grot1=grot1/max(abs(grot1(:)));
for i=1:size(grot1,1)
  for j=1:size(grot1,1)
    if clustersALL(grotDD(i)) == clustersALL(grotDD(j)) % why doing this?
      grot1(i,j)=grot1(i,j)+1;
    end
  end
end


%[hier,linkages]=nets_hierarchy(grot*300,Znet1, [1:200],SUMPICS,1.5); 

% To get the "grid" of edges:
%grot=ZnetMOD;
% nets_edgepics(ts,SUMPICS,Mnet2,grot,200,2); 
% set(gcf,'PaperPositionMode','auto','Position',[10 10 1800 900]); 
%print('-dpng',sprintf('%s/edgemod.png','/home/fs0/steve'));

[hier,linkages]=nets_hierarchy(grot1,grot*10,grotDD,SUMPICS,0.75); 
  set(gcf,'PaperPositionMode','auto','Position',[10 10 2800 2000]);   %print('-dpng',sprintf('%s/edgemodhier.png','/home/fs0/steve'));
clusters=cluster(linkages,'maxclust',4)';

% system(sprintf('/bin/rm -rf %s/netjs',NM));
% system(sprintf('cp -r %s/netjs %s',FMRIB,NM));
NP=sprintf('%s/netjs/data/dataset1',NM);
save(sprintf('%s/Znet3.txt',NP),'grot','-ascii');
grot=grot.*sign(Mnet2(grotDD,grotDD));
save(sprintf('%s/Znet4.txt',NP),'grot','-ascii');
grot=Mnet1(grotDD,grotDD);   save(sprintf('%s/Znet1.txt',NP),'grot','-ascii');
grot=Mnet2(grotDD,grotDD);   save(sprintf('%s/Znet2.txt',NP),'grot','-ascii');
save(sprintf('%s/hier.txt',NP),'hier','-ascii');
save(sprintf('%s/linkages.txt',NP),'linkages','-ascii');
save(sprintf('%s/clusters.txt',NP),'clusters','-ascii');
system(sprintf('/bin/mkdir %s/melodic_IC_sum.sum',NP));
for i=1:length(grotDD)
  system(sprintf('/bin/cp %s/%.4d.png %s/melodic_IC_sum.sum/%.4d.png',SUMPICS,grotDD(i)-1,NP,i-1));
end
! cp ~/main.js ~/netvis.js ~/www/HCP_GigaTrawl/netjs/js

% run following in terminal
% open /Applications/Google\ Chrome.app --args --allow-file-access-from-files

%%
% pos_map = ZnetMOD>0;
% neg_map = ZnetMOD<0;
% 
% for i = 1:200
%     Sum_pos(i) = sum(ZnetMOD(i,pos_map(i,:)));
%     Sum_neg(i) = sum(ZnetMOD(i,neg_map(i,:)));
% end
% 
% for i = 1:200
%     mean_pos(i) = mean(ZnetMOD(i,pos_map(i,:)));
%     mean_neg(i) = mean(ZnetMOD(i,neg_map(i,:)));
% end
% 
% figure
% yyaxis left
% plot(Sum_pos)
% yyaxis right
% plot(mean_pos)
% legend('sum','mean')
% title('Sum vs mean positive map')
% 
% figure
% yyaxis left
% plot(Sum_neg)
% yyaxis right
% plot(mean_neg)
% legend('sum','mean')
% title('Sum vs mean negative map')


% Fisher's z transformation and multiplied by sign of mean map
%ZnetMOD_z=0.5*log((1+ZnetMOD)./(1-ZnetMOD))./std(ZnetMOD,'omitnan'); %an alternative of the next line
ZnetMOD_z=(1/mean(std(ZnetMOD,'omitnan')))*0.5*log((1+ZnetMOD)./(1-ZnetMOD)); % multiply by a constant is to preserve the variance in std between nodes. 22 is obtained by 1/mean(std(all nodes)).
grot = ZnetMOD_z;
grot=grot.*sign(Mnet2);
grot = sort(grot);

%%% sum up top 50 positive and negative weights to form pos and neg maps
Sum_neg = sum(-grot(1:20,:));
Sum_pos = sum(grot(181:200,:));

figure
plot(Sum_pos)
hold on 
plot(Sum_neg)
legend('Positive map','Negative map')
title('Normalised by node constant')

figure 
scatter(Sum_pos, Sum_neg)
hold on
refline([1,0])
xlabel('Positive loadings')
ylabel('Negative loadings')
% pos_map = grot>0;
% neg_map = grot<0;
% 
% for i = 1:200
%     Sum_pos(i) = sum(grot(i,pos_map(i,:)));
%     Sum_neg(i) = sum(grot(i,neg_map(i,:)));
% end

% for i = 1:200
%     mean_pos(i) = mean(grot(i,pos_map(i,:)));
%     mean_neg(i) = mean(grot(i,neg_map(i,:)));
% end

figure
yyaxis left
plot(Sum_pos)
yyaxis right
plot(mean_pos)
legend('sum','mean')
title('Sum vs mean positive map')

figure
yyaxis left
plot(Sum_neg)
yyaxis right
plot(mean_neg)
legend('sum','mean')
title('Sum vs mean negative map')


writetable(array2table(Sum_pos'), '../BrainVis/CCA1_pos_unmapped.txt', 'WriteVariableNames',false);
writetable(array2table(Sum_neg'), '../BrainVis/CCA1_neg_unmapped.txt', 'WriteVariableNames',false);
%%
WBC = '/Applications/workbench/bin_macosx64/wb_command';
BO=ciftiopen('/Users/jessie_liu/Desktop/HCP1200/rawdata/groupICA/melodic_IC.dscalar.nii',WBC); 
GM=BO.cdata;
GM=GM.*repmat(sign(max(GM)+min(GM)),size(GM,1),1)./repmat(max(abs(GM)),size(GM,1),1);  % make all individual group maps have a positive peak, and of peak height=1

BO.cdata=log(mean(-GM * sum(grot(1:20,:))',2)    ./ mean(GM,2)); 
ciftisavereset(BO,'SDRCCA1_NEG20_unmapped.dscalar.nii',WBC);    
%[prctile(BO.cdata,80) corr(BO.cdata(1:59000,1),clus2(1:59000,1))]
BO.cdata=log(mean( GM * sum(grot(181:200,:))',2) ./ mean(GM,2)); 
ciftisavereset(BO,'SDRCCA1_POS20_unmapped.dscalar.nii',WBC);    

[prctile(BO.cdata,80) corr(BO.cdata(1:59000,1),clus2(1:59000,1))]
%BO.cdata=log(mean( GM * sum(grot)',2)            ./ mean(GM,2)); ciftisave(BO,'~/sumCCAnodesPOSNEG.dscalar.nii',WBC);