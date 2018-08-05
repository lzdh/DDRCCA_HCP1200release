clear
clc
behaviour = importdata('../Data/SM1003x387.mat');
headers = importdata('../Data/headers_SM387.txt');


%%% New grouping
confounds = [2:15]; %1
Demographics = [16:22];%2
PhyHeal = [23:36];%3
FemHeal = [37:45];%4
FamHist = [46:63];%5
Psychiatry = [64:107];%6
Sensory = [108:125];%7
Drug = [126:143];%8
Alcohol = [144:177];%9
Tobacco = [178:213];%10
Alertness = [214:239];%11
Cognition = [240:291];%12
Emotion = [292:315];%13
Motor = [316:322];%14
Personality = [323:387];%15


% remove ill-conditioned variables
sub_area = {confounds, Demographics,PhyHeal,FemHeal,FamHist,Psychiatry, ...
            Sensory,Drug,Alcohol,Tobacco,Alertness,Cognition,Emotion,Motor,Personality};
headers_new = headers{1};
behaviour_new = behaviour(:,1);
indi = 0;
for d = 1:15
    badvars=[]; %exclude bad variables in the sub-area
    SM_sub=behaviour(:,sub_area{d});
    header_d=headers(sub_area{d});
    
    for i=1:size(SM_sub,2)
      Y=SM_sub(:,i); grotKEEP=~isnan(Y);  
      %change the creterion when # of subjects changes
      if d==4
          if (sum(grotKEEP)>=250) & (std(Y(grotKEEP))>0)& max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95 % & (grot<100);
          % the 3rd thing above is:  is the size of the largest equal-values-group too large?
            i=i;
          else
              [i sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))];
              badvars=[badvars i];
          end
      else
          if (sum(grotKEEP)>500) & (std(Y(grotKEEP))>0)& (max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))<0.95)
            % the 3rd thing above is:  is the size of the largest equal-values-group too large?
            i=i; % do nothing
          else
            [i sum(grotKEEP)  max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP))]
            badvars=[badvars i];
          end 
      end
    end
    varskeep = setdiff([1:size(sub_area{d},2)],badvars);
    behaviour_new = [behaviour_new, SM_sub(:,varskeep)];
    sub_area{d} = [(size(behaviour_new,2)-length(varskeep)+1):size(behaviour_new,2)];% update sub-domain indices
    headers_new = [headers_new; header_d(varskeep,:)];
end

behaviour = behaviour_new;
headers = headers_new;


% calculate pairwise correlation

varsd=palm_inormal(behaviour);
% varsdCOVt=zeros(size(varsd,2));
% for i=1:size(varsd,2)
%   for j=1:size(varsd,2)
%     grot=varsd(:,[i j]); 
%     grot=corrcoef(grot(sum(isnan(grot),2)==0,:)); 
%     varsdCOVt(i,j)=grot(1,2);
%   end
% end
% varsdCOVt2=nearestSPD(varsdCOVt); % scatter(varsdCOVt(:),varsdCOVt2(:));
% % Ensure it's a correlation matrix!
% varsdCOVt2 = diag(sqrt(diag(varsdCOVt2)).^(-1/2))*varsdCOVt2*diag(sqrt(diag(varsdCOVt2)).^(-1/2));
% scatter((varsdCOVt(:)+varsdCOVt2(:))/2,varsdCOVt(:)-varsdCOVt2(:))
% % Avoid embaressing 'nearly 1' diagonal entries
% varsdCOVt2(diag(diag(ones(size(varsdCOVt2))))>0)=1;
% 
% %imagesc(varsdCOVt2);axis image; colorbar
figure
imagesc(corr(behaviour, 'rows', 'pairwise'));axis image; colorbar
title('Correlation among 318 unflipped behavioral and demongraphic measures')


%% this doesnt really work.... 
dd=mean(varsdCOVt2.*(abs(varsdCOVt2)>0.0));  
headers_flipped = headers;
varsdCOVt2f = varsdCOVt;
vars_flip = behaviour;

for i = find(dd<0)
  varsdCOVt2f(i,:)=-varsdCOVt2f(i,:);
  varsdCOVt2f(:,i)=-varsdCOVt2f(:,i);
  vars_flip(:,i)=-vars_flip(:,i);
  headers_flipped{i} = ['-' headers_flipped{i}];
end

figure(1);imagesc(varsdCOVt,[-.5 .75]);axis image; colorbar
figure(2);imagesc(varsdCOVt2f,[-.5 .75]);axis image; colorbar

%% try: set a standard variable and all correlate with this one
stdd = varsd(:,18); %income
headers_flipped = headers;
vars_flip = behaviour;
for i =1:size(varsd,2)
    vars_corr(i) = corr(varsd(:,i),stdd, 'rows', 'pairwise');
end

vars_corr(1,[sub_area{9},sub_area{11},82,99,180,202,211]) = -0.2; %flip
vars_corr(1,[1,sub_area{1},sub_area{2}, sub_area{15}, 25:30,126:128,130,132:134,146,156, 85, 86]) = 0.2; %dont flip

for i =1:size(varsd,2)
    if vars_corr(i)<0
        vars_flip(:,i) = -vars_flip(:,i);
        headers_flipped{i} = ['-' headers_flipped{i}];  
    end
end


%Plot
%%%Define centriod for each domain
for i =1:15
    centriod(i) = sub_area{i}(1);
end   

figure
imagesc(corr(vars_flip, 'rows', 'pairwise'));axis image; colorbar
title('Correlation among 318 flipped behavioral and demongraphic measures')
ax = gca;
ax.XTick = centriod;
ax.XTickLabel = {'confounds', 'Demographics','PhyHeal','FemHeal','FamHist','Psychiatry', ...
            'Sensory','Drug','Alcohol','Tobacco','Alertness','Cognition','Emotion','Motor','Personality'};
ax.XTickLabelRotation = 45;
ax.YTick = centriod;
ax.YTickLabel = {'confounds', 'Demographics','PhyHeal','FemHeal','FamHist','Psychiatry', ...
            'Sensory','Drug','Alcohol','Tobacco','Alertness','Cognition','Emotion','Motor','Personality'};

writetable(cell2table(headers_flipped), 'header318_flipped.txt', 'WriteVariableNames',false);

%% without subquestions
confounds = [2:15]; %1
Demographics = [16:22];%2
PhyHeal = [23:30];%3
FemHeal = [31:35];%4
FamHist = [36:40];%5
Psychiatry = [41:84];%6
Sensory = [85:96];%7
Drug = [97:107];%8
Alcohol = [108:135];%9
Tobacco = [136:145];%10
%Alertness = [146:171];%11
Alertness = [146:154];%11
%Cognition = [172:223];%12
Cognition = [172:215];%12
Emotion = [224:246];%13
Motor = [247:253];%14
%Personality = [254:318];%15
Personality = [254:258];%15

sub_area = {confounds, Demographics,PhyHeal,FemHeal,FamHist,Psychiatry, ...
            Sensory,Drug,Alcohol,Tobacco,Alertness,Cognition,Emotion,Motor,Personality};

vars_flip_short = vars_flip(:,[1,confounds, Demographics,PhyHeal,FemHeal,FamHist,Psychiatry, ...
            Sensory,Drug,Alcohol,Tobacco,Alertness,Cognition,Emotion,Motor,Personality]);

headers_flipped_short = headers_flipped([1,confounds, Demographics,PhyHeal,FemHeal,FamHist,Psychiatry, ...
            Sensory,Drug,Alcohol,Tobacco,Alertness,Cognition,Emotion,Motor,Personality]);

centriod = [2,16,23,31,36,41,85,97,108,136,146,155,199,222,229];      
        
figure
imagesc(corr(vars_flip_short, 'rows', 'pairwise'));axis image; colorbar
title('Correlation among 318 flipped behavioral/demongraphic measures grouped by domains')
ax = gca;
ax.XTick = centriod;
ax.XTickLabel = {'confounds', 'Demographics','PhyHeal','FemHeal','FamHist','Psychiatry', ...
            'Sensory','Drug','Alcohol','Tobacco','Alertness','Cognition','Emotion','Motor','Personality'};
ax.XTickLabelRotation = 90;
ax.YTick = centriod;
ax.YTickLabel = {'confounds', 'Demographics','PhyHeal','FemHeal','FamHist','Psychiatry', ...
            'Sensory','Drug','Alcohol','Tobacco','Alertness','Cognition','Emotion','Motor','Personality'};

writetable(cell2table(headers_flipped_short), 'header233_flipped.txt', 'WriteVariableNames',false);


 