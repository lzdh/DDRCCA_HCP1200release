clear
%clc
SM = importdata('../Data/Flipped1003x233.mat');% to include subquestions load 'Flipped1003x318.mat'
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

%%% setup confounds matrix
conf=palm_inormal([SM(:,confounds(1:12)), SM(:,confounds(13:14)).^(1/3)]); % Gaussianise.  479 480 481 race variables
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf, conf(:,5).^2, conf(:,10:12).^2]);  % add on squared terms and renormalise

sub_area={Demographics,PhyHeal,FemHeal, FamHist, Psychiatry, Sensory, Drug,Alcohol, Tobacco, Alertness, Cognition, Emotion, ...
             Motor, Personality};
         
[C0,IA]=unique(Fam_info.Family_ID);
PC = [];

for d=1:length(sub_area)
    fprintf('\n Estimating dimension for sub-area %d', d);
    SM_sub = SM(:,sub_area{d});
    
    varsd=palm_inormal(SM_sub); % Gaussianise
    for i=1:size(varsd,2) % deconfound ignoring missing data
      grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
    end
    
    %%% Estimating dimension (minimising  PRESS) for the sub-domain using
    %%% 5-fold CV for 10 repetitions
    
    for k=1:20
        fprintf('\n CV repetition %d', k);
        
        error2=0;
        error1=0;
        Error1_sim = [];
        Error2_sim = [];
        
        K=5;
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
            
%             grot=varsd_in; grotI=double(~isnan(grot)); grot(isnan(grot))=0;
%             varsdCOV = (grot'*grot) ./ (grotI'*grotI);
%             varsdCOV=nearestSPD(varsdCOV);
%             [vv1,dd1]=eigs(varsdCOV,size(varsd_in,2));
             
            COV=zeros(size(varsd_in,2));
            for i=1:size(varsd_in,2) % estimate "pairwise" covariance, ignoring missing data
              for j=1:size(varsd_in,2)
                grot=varsd_in(:,[i j]); grot=cov(grot(sum(isnan(grot),2)==0,:)); COV(i,j)=grot(1,2);
              end
            end
            varsdCOV=nearestSPD(COV);
            [vv1,dd1]=eigs(varsdCOV,size(varsd_in,2));
            
            
            if dd1(1,1)<dd1(2,2)  vv1=fliplr(vv1); end
            varsd_out(isnan(varsd_out))=0; % impute missing data as zeros
            
            err2=[];
            for j=1:min(size(vv1,2),20)
                P = vv1(:,1:j)*vv1(:,1:j)';
                err1 = varsd_out * (eye(size(P)) - P); % err1 is needed for sub-domain report
                
                for i=1:size(varsd_out,2)
                    proj = varsd_out(:,[1:i-1 i+1:end])*pinv(vv1([1:i-1 i+1:end],1:j))'*vv1(:,1:j)'; 
                    err2(:,i) = varsd_out(:,i) - proj(:,i); %Pseudoinverse
                end
                error1(index2,j) = sum(err1,2).^2;
                error2(index2,j) = sum(err2,2).^2;
            end
        end
        Error1_sim(k,:,:) = error1;
        Error2_sim(k,:,:) = error2;
        temp = sum(error2);
        N_sim(k,d) = find(temp==(min(temp)));
    end
    Error1_ave = reshape(mean(Error1_sim,1), [size(error1,1),size(error1,2)]);
    Error2_ave = reshape(mean(Error2_sim,1), [size(error1,1),size(error1,2)]);
    
    Error1{d} = sum(Error1_ave); %Error1 is needed for sub-domain report
    Error2{d} = sum(Error2_ave); 
    
    N_pca(d)=find(Error2{d}==(min(Error2{d})));
    
    % calculate principal components
    varsd_t=varsd;
    varsd_t(isnan(varsd_t))=0;
    [u,s,v]=svd(varsd_t);
    PC = [PC, u*s(:,1:N_pca(d))];
    
end

%print total number of variables in each sub-area and the estimated dimension
cellfun(@length, sub_area)
N_pca

%% Brain measures: prepare brain measrues for SDR
NETmat = importdata('../rawdata/netmats2.txt');

%%% deconfound netmat
NET1=nets_demean(NETmat);  NET1=NET1/std(NET1(:)); % no norm
NETd=nets_demean(NET1-conf*(pinv(conf)*NET1));

NETd_nodes = reshape(NETd, [1003,200,200]);

for i=1:200
    Node=[];
    for n=1:size(NETmat,1)
        temp=NETd_nodes(n,i,:);
        Node=[Node; temp];
    end
    NETd_node{i}=reshape(Node,[size(NETmat,1) 200]);
end

%% Two-way CV on brain measures
Fam_info=readtable('../Data/Fam_infoBMsubs.csv');
[C0,IA]=unique(Fam_info.Family_ID);
PC_bm = [];
SM = importdata('../Data/Flipped1003x233.mat');% to include subquestions load 'Flipped1003x318.mat'
K=5; % K-fold CV
load('../Data/NETd_node.mat');

%%% setup confounds matrix
confounds = [2:15]; %0
conf=palm_inormal([SM(:,confounds(1:12)), SM(:,confounds(13:14)).^(1/3)]); % Gaussianise.  479 480 481 race variables
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf, conf(:,5).^2, conf(:,10:12).^2]);  % add on squared terms and renormalise


for d=1:200
    fprintf('\n Estimating dimension for node %d', d);
    
    Node=NETd_node{d};

    %%% normalise each node
    Noded=nets_demean(Node);  Noded=Noded/std(Noded(:)); % no norm
    Noded=nets_demean(Noded-conf*(pinv(conf)*Noded));   % deconfound and demean

    
    for k=1:5
        fprintf('\n CV repetition %d', k);
        
        error1=0;
        Error2_sim = [];
        
        groups = crossvalind('Kfold',length(C0),K); %K-fold
        Fam_info.group_ID=zeros(size(SM,1),1);
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
            Noded_out=Noded(index2',:);

            for i=1:size(Train,1)
                index1(i)=find(Train(i)==SM(:,1));
            end
            Noded_in=Noded(index1',:);
               
            [~,dd1,vv1]=svd(Noded_in);
            if dd1(1,1)<dd1(2,2)  vv1=fliplr(vv1); end
               
            err2=[];
            for j=1:60              
                for i=1:size(Noded_out,2)
                    proj = Noded_out(:,[1:i-1 i+1:end])*pinv(vv1([1:i-1 i+1:end],1:j))'*vv1(:,1:j)'; 
                    err2(:,i) = Noded_out(:,i) - proj(:,i); %Pseudoinverse
                end
                error2(index2,j) = sum(err2,2).^2;
            end
        end
        Error2_sim(k,:,:) = error2;
    end
    Error2_ave = reshape(mean(Error2_sim,1), [size(error2,1),size(error2,2)]);    
    Error2{d} = sum(Error2_ave); 
    N_pca(d)=find(Error2{d}==(min(Error2{d})));
    
    % calculate principal components
%     [u,s,v]=svd(Noded);
%     PC_bm = [PC_bm, u*s(:,1:N_pca(d))];
     
end

%% Brain measures buster (cluster) results
clear
clc

N=[];
for k=1:200
  myfilename = sprintf('/Users/zliu/Desktop/HCP1200/Result/buster_results/batch2/SDRNET_node%d.mat', k);
  load(myfilename);
  Error{k} = Error2;
  N=[N N_pca];
end
 
load('../Data/NETd_node.mat')

PC_bm=[];
for i=1:200
   Netmat=NETd_node{i}; 
  [uu,ss,vv]=nets_svds(Netmat,N(i)); % SVD reduction
   pc=uu*ss;
  PC_bm=[PC_bm pc];
end

for i =1:200
    plot(Error{i})
    hold on
end
