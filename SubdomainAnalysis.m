clear
clc
SM = importdata('../Data/Flipped1003x233.mat');
headers = importdata('../Data/header233_flipped.txt'); 
Fam_info=readtable('../Data/Fam_infoBMsubs.csv');

% Define sub-domains
confounds = [2:15]; %0
Demographics = [16:22];%1
PhyHeal = [23:30];%2
FemHeal = [31:35];%3
FamHist = [36:40];%4
Psychiatry = [41:84];%5
Sensory = [85:96];%6
Drug = [97:107];%7
Alcohol = [108:135];%8
Tobacco = [136:145];%9
%Alertness = [146:171];%10
Alertness = [146:154];%10
%Cognition = [172:223];%11
Cognition = [155:198];%11
Emotion = [199:221];%12
Motor = [222:228];%13
%Personality = [254:318];%14
Personality = [229:233];%14

sub_dom = SM(:,Personality);
sub_header = headers(Personality);

%%% setup confounds matrix
conf=palm_inormal([SM(:,confounds(1:12)), SM(:,confounds(13:14)).^(1/3)]); % Gaussianise.  479 480 481 race variables
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf, conf(:,5).^2, conf(:,10:12).^2]);  % add on squared terms and renormalise

varsd=palm_inormal(sub_dom); % Gaussianise
for i=1:size(varsd,2) % deconfound ignoring missing data
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end

%------------------------------------------------------------------------%
% Eigen-spectrum
%------------------------------------------------------------------------%

COV_pca=zeros(size(varsd,2));
for i=1:size(varsd,2) % estimate "pairwise" covariance, ignoring missing data
  for j=1:size(varsd,2)
    grot=varsd(:,[i j]); grot=cov(grot(sum(isnan(grot),2)==0,:)); COV_pca(i,j)=grot(1,2);
  end
end
varsdCOV=nearestSPD(COV_pca);
[uu,dd]=eigs(varsdCOV,size(varsd,2)); 

eigen_spec=sort(diag(dd),'descend')';

%------------------------------------------------------------------------%
% Null eigen-spectrum
%------------------------------------------------------------------------%

for j = 1:size(varsd,2)
    p = randperm(size(varsd,1));
    varsd_n(:,j) = varsd(p,j);
end

% varsd_n(isnan(varsd_n))=0;
% [u_n,d_n,v_n]=svd(varsd_n);
% d_n = d_n.^2;

COV_n=cov(varsd_n,'partialrows');
varsdCOV2=nearestSPD(COV_n);
[u_n,d_n]=eigs(varsdCOV2,size(varsd_n,2)); 

null_eigen_spec=sort(diag(d_n),'descend')';

%------------------------------------------------------------------------%
% CV principal components
%------------------------------------------------------------------------%
[C0,IA]=unique(Fam_info.Family_ID);        
K=5;

%for k = 1:10
groups = crossvalind('Kfold',length(C0),K); %K-fold
Fam_info.group_ID=zeros(size(varsd,1),1);
for i=1:size(Fam_info,1)
    index=find(strcmp(Fam_info.Family_ID{i},C0));
    Fam_info.group_ID(i)=groups(index);
end          

for n=1:K

    Test=Fam_info(Fam_info.group_ID==n,1);
    Test=table2array(Test);
    Train=setdiff(Fam_info.Subject,Test);

    index1=0;
    index2=0;
    for i=1:size(Test,1)
        index2(i)=find(Test(i)==SM(:,1));
    end
    varsd_out=varsd(index2',:);

    for i=1:size(Train,1)
        index1(i)=find(Train(i)==SM(:,1));
    end
    varsd_in=varsd(index1',:);

    varsd1=palm_inormal(varsd_in); %Train set
    varsd2=palm_inormal(varsd_out); % Test set
    
    varsd1(isnan(varsd1))=0;
    varsd2(isnan(varsd2))=0;
    
    varsdCOV=varsd1'*varsd1;
    [vv1,dd1]=eigs(varsdCOV,size(varsd1,2));  % SVD (eigs actually)
    
    if dd1(1,1)<dd1(end,end) vv1=fliplr(vv1);end
       
    varsdCOV2=varsd2'*varsd2;
    var2=trace(varsdCOV2);

    for j=1:(size(varsdCOV,1))
        
        Test_pca=varsd2*vv1(:,1:j);
        Test_cov=Test_pca'*Test_pca;
    
        r(n,j)=trace(Test_cov)/var2;
    end

    [vv_r(:,:,n)]=rotatefactors(vv1(:,1:3)); % take the first n PCs and apply the varimax
    [vv_obl(:,:,n)]=rotatefactors(vv1(:,1:3), 'Method', 'promax');
    %    %inv(T'*T)
    v1(:,:,n)=vv1(:,1:3);

    eigen_spec_CV(n,:)=sort(diag(dd1),'descend')';

end

rbar=mean(r);

for j=1:size(v1,2)
    for k=1:K
        Corr=corr(v1(:,j,1),v1(:,j,k));
        if Corr<0 v1(:,j,k)=-v1(:,j,k); end
            
        Corr=corr(vv_r(:,j,1),vv_r(:,j,k));
        if Corr<0 vv_r(:,j,k)=-vv_r(:,j,k); end
            
        Corr=corr(vv_obl(:,j,1),vv_obl(:,j,k));
        if Corr<0 vv_obl(:,j,k)=-vv_obl(:,j,k); end
    end
end

for j=1:size(vv_r,2)
    for k=1:K
        w(k,:)=v1(:,j,k);
        w_r(k,:)=vv_r(:,j,k);
        w_olb(k,:)=vv_obl(:,j,k);
    end
    weight(:,j)=mean(w);
    weight_r(:,j)=mean(w_r);
    weight_obl(:,j)=mean(w_olb);
end


%%% Flip the weights to make sure mean weights are positive
for j=1:size(vv_r,2)
    tmp = mean(weight(:,j));
    tmp_r = mean(weight_r(:,j));
    tmp_obl = mean(weight_obl(:,j));
    if tmp<0
        weight(:,j) = -weight(:,j);
    end
   
    if tmp_r<0
        weight_r(:,j) = -weight_r(:,j);
    end
    if tmp_obl<0
        weight_obl(:,j) = -weight_obl(:,j);
    end
   
end

