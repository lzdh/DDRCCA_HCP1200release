clear
clc

SMvars = importdata('../Data/Flipped1003x233.mat');
headers = importdata('../Data/header233_flipped.txt'); 
NET = importdata('../rawdata/netmats2.txt');
% N_SM = importdata('../Data/N_pca_SM.mat');
% N_BM = importdata('../Data/N_pca_BM.mat');
load('../Data/NETd_node.mat');
Fam_info=readtable('../Data/Fam_infoBMsubs.csv');
%SM_sdr = importdata('../Data/SDR_SM62.mat');
%BM_sdr = importdata('../Data/SDR_BM6867.mat');
EB=hcp2blocks('../rawdata/restricted.txt', [ ], false, SMvars(:,1));


%%% Generate SDR variables and do factor rotation
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

sub_area={Demographics,PhyHeal,FemHeal, FamHist, Psychiatry, Sensory, Drug,Alcohol, Tobacco, Alertness, Cognition, Emotion, ...
             Motor, Personality};

%%% setup confounds matrix
% conf=palm_inormal([SMvars(:,confounds(1:12)), SMvars(:,confounds(13:14)).^(1/3)]); 
% conf(isnan(conf))=0;  % impute missing data as zeros
% conf=nets_normalise([conf, conf(:,5).^2, conf(:,10:12).^2]);  % add on squared terms and renormalise
% 
% sub_area={Demographics,PhyHeal,FemHeal, FamHist, Psychiatry, Sensory, Drug,Alcohol, Tobacco, Alertness, Cognition, Emotion, ...
%              Motor, Personality};
% 
% %%% deconfound and normalise original SM and BM
% varsd_train=palm_inormal(SMvars);
% for i=1:size(varsd_train,2)
%   grot=(isnan(varsd_train(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd_train(grot,i)=nets_normalise(varsd_train(grot,i)-grotconf*(pinv(grotconf)*varsd_train(grot,i)));
% end
% 
% NET1=nets_normalise(NET); % no norm
% NETd=nets_normalise(NET1-conf*(pinv(conf)*NET1));   % deconfound and demean
Nperm=10000;
[C0,IA]=unique(Fam_info.Family_ID);
K=5;
%for k=1:20
    
    %fprintf('\n Simulation # %d',k);
    

groups = crossvalind('Kfold',length(C0),K); %K-fold
Fam_info.group_ID=zeros(size(SMvars,1),1);
for i=1:size(Fam_info,1)
    index=find(strcmp(Fam_info.Family_ID{i},C0));
    Fam_info.group_ID(i)=groups(index);
end          

for n=1:K
    fprintf('\n Fold %d', n);
    Test=Fam_info(Fam_info.group_ID==n,1);
    Test=table2array(Test);
    Train=setdiff(Fam_info.Subject,Test);

    index_train=0;
    index_test=0;
    for i=1:size(Test,1)
        index_test(i)=find(Test(i)==SMvars(:,1));
    end
    vars_test=SMvars(index_test,:);
    net_test = NET(index_test,:);

    for i=1:size(Train,1)
        index_train(i)=find(Train(i)==SMvars(:,1));
    end
    vars_train=SMvars(index_train,:);
    net_train = NET(index_train,:);

    %%% construct confounds for training and test sets
    conf_train=palm_inormal([vars_train(:,confounds(1:12)), vars_train(:,confounds(13:14)).^(1/3)]); 
    conf_train(isnan(conf_train))=0;  % impute missing data as zeros
    conf_train=nets_normalise([conf_train, conf_train(:,5).^2, conf_train(:,10:12).^2]);  % add on squared terms and renormalise

    conf_test=palm_inormal([vars_test(:,confounds(1:12)), vars_test(:,confounds(13:14)).^(1/3)]); 
    conf_test(isnan(conf_test))=0;  % impute missing data as zeros
    conf_test=nets_normalise([conf_test, conf_test(:,5).^2, conf_test(:,10:12).^2]);  % add on squared terms and renormalise

    %%% normalise and de-confound SM
    varsd_train=palm_inormal(vars_train); % Gaussianise
    for i=1:size(varsd_train,2) % deconfound ignoring missing data
      grot=(isnan(varsd_train(:,i))==0); grotconf=nets_demean(conf_train(grot,:)); varsd_train(grot,i)=normalise(varsd_train(grot,i)-grotconf*(pinv(grotconf)*varsd_train(grot,i)));
    end

    varsd_test=palm_inormal(vars_test); % Gaussianise
    for i=1:size(varsd_test,2) % deconfound ignoring missing data
      grot=(isnan(varsd_test(:,i))==0); grotconf=nets_demean(conf_test(grot,:)); varsd_test(grot,i)=normalise(varsd_test(grot,i)-grotconf*(pinv(grotconf)*varsd_test(grot,i)));
    end

    %%% normalise and de-confound BM
    netd_train = nets_demean(net_train); netd_train = netd_train/std(netd_train(:));
    netd_train = nets_demean(netd_train-conf_train*(pinv(conf_train)*netd_train));

    netd_test = nets_demean(net_test); netd_test = netd_test/std(netd_test(:));
    netd_test = nets_demean(netd_test-conf_test*(pinv(conf_test)*netd_test));


    %%% for training set, apply two-way CV to estimate dims for all for
    %%% both SM and BM sub-domains (this requires another 5-fold CV)
    Fam_info_sub = Fam_info(index_train,:);
    [C1,~]=unique(Fam_info_sub.Family_ID);

    SDRgroups = crossvalind('Kfold',length(C1),K); %K-fold
    Fam_info_sub.group_ID=zeros(size(vars_train,1),1);
    for i=1:size(Fam_info_sub,1)
        index=find(strcmp(Fam_info_sub.Family_ID{i},C1));
        Fam_info_sub.group_ID(i)=SDRgroups(index);
    end 

    Error_sm=[];
    Error_bm=[];
    for m=1:K
        
        held_out=Fam_info_sub(Fam_info_sub.group_ID==m,1);
        held_out=table2array(held_out);
        held_in=setdiff(Fam_info_sub.Subject,held_out);

        index1=0;
        index2=0;
        for i=1:size(held_out,1)
            index2(i)=find(held_out(i)==vars_train(:,1))';
        end

        for i=1:size(held_in,1)
            index1(i)=find(held_in(i)==vars_train(:,1))';
        end
        vars_in=vars_train(index1,:);
        vars_out=vars_train(index2,:);

        %%% construct confounds for held-in and held-out sets
        conf_in=palm_inormal([vars_in(:,confounds(1:12)), vars_in(:,confounds(13:14)).^(1/3)]); 
        conf_in(isnan(conf_in))=0;  % impute missing data as zeros
        conf_in=nets_normalise([conf_in, conf_in(:,5).^2, conf_in(:,10:12).^2]);  % add on squared terms and renormalise

        conf_out=palm_inormal([vars_out(:,confounds(1:12)), vars_out(:,confounds(13:14)).^(1/3)]); 
        conf_out(isnan(conf_out))=0;  % impute missing data as zeros
        conf_out=nets_normalise([conf_out, conf_out(:,5).^2, conf_out(:,10:12).^2]);  % add on squared terms and renormalise

        %%% normalise and de-confound SM
        varsd_in=palm_inormal(vars_in); % Gaussianise
        for i=1:size(varsd_in,2) % deconfound ignoring missing data
          grot=(isnan(varsd_in(:,i))==0); grotconf=nets_demean(conf_in(grot,:)); varsd_in(grot,i)=normalise(varsd_in(grot,i)-grotconf*(pinv(grotconf)*varsd_in(grot,i)));
        end

        varsd_out=palm_inormal(vars_out); % Gaussianise
        for i=1:size(varsd_out,2) % deconfound ignoring missing data
          grot=(isnan(varsd_out(:,i))==0); grotconf=nets_demean(conf_out(grot,:)); varsd_out(grot,i)=normalise(varsd_out(grot,i)-grotconf*(pinv(grotconf)*varsd_out(grot,i)));
        end

        %%% SDR on SM
        for d=1:length(sub_area)
            fprintf('\n Estimating dimension for sub-area %d', d);

            varsd_in_d = varsd_in(:,sub_area{d});
            varsd_out_d  = varsd_out(:,sub_area{d});

            COV_in = cov(varsd_in_d,'partialrows');
            varsdCOV=nearestSPD(COV_in);
            [vv1,dd1]=eigs(varsdCOV,size(varsd_in_d,2));

            if dd1(1,1)<dd1(2,2)  vv1=fliplr(vv1); end
            varsd_out_d(isnan(varsd_out_d))=0; % impute missing data as zeros

            err2=[];
            error2=[];
            for j=1:min(size(vv1,2),20)
                for i=1:size(varsd_out_d,2)
                    proj = varsd_out_d(:,[1:i-1 i+1:end])*pinv(vv1([1:i-1 i+1:end],1:j))'*vv1(:,1:j)'; 
                    err2(:,i) = varsd_out_d(:,i) - proj(:,i); %Pseudoinverse
                end
                error2(:,j) = sum(err2,2).^2;
            end
            Error_sm{d}(index2,:) = error2;
        end
        
        
        %%% SDR on BM
        for d=1:200
            net_in = NETd_node{d}(index_train,:);
            fprintf('\n Estimating dimension for node %d', d);

            net_in_d = net_in(index1,:);
            net_out_d = net_in(index2,:);

            %de-mean normailise by node std (weight node variance) and de-confound held-in and held-out sets
            net_in_d = nets_demean(net_in_d); net_in_d = net_in_d/std(net_in_d(:));
            net_out_d = nets_demean(net_out_d); net_out_d = net_out_d/std(net_out_d(:));

            net_in_d=nets_demean(net_in_d-conf_in*(pinv(conf_in)*net_in_d));
            net_out_d=nets_demean(net_out_d-conf_out*(pinv(conf_out)*net_out_d));

            [~,dd2,vv2]=svd(net_in_d);
            if dd2(1,1)<dd2(2,2)  vv2=fliplr(vv2); end

            err=[];
            error=[];
            for j=1:40              
                for i=1:size(net_out_d,2)
                    proj = net_out_d(:,[1:i-1 i+1:end])*pinv(vv2([1:i-1 i+1:end],1:j))'*vv2(:,1:j)'; 
                    err(:,i) = net_out_d(:,i) - proj(:,i); %Pseudoinverse
                end
                error(:,j) = sum(err,2).^2;
            end
            Error_bm{d}(index2,:) = error;
        end    
    end

    %%%construct latent factors for SM and BM to feed into CCA
    rPC_SM_train = [];
    rPC_SM_test = [];
    for d = 1:length(sub_area)
        temp = sum(Error_sm{d});
        N_sm(n,d) = find(temp==(min(temp))); 

        varsd_train(isnan(varsd_train))=0;
        varsd_test(isnan(varsd_test))=0;

        COV = cov(varsd_train(:,sub_area{d}), 'partialrows');
        [v,~] = eigs(COV, N_sm(n,d));
        
        if N_sm(n,d) ~= 1
%             [v_r,T] = rotatefactors(v);
%             rPC_SM_train = [rPC_SM_train, varsd_train*v_r];
%             rPC_SM_test = [rPC_SM_test, varsd_test*v_r];
           
            rPC_SM_train = [rPC_SM_train, varsd_train(:,sub_area{d})*v];
            rPC_SM_test = [rPC_SM_test, varsd_test(:,sub_area{d})*v];
        else
            rPC_SM_train = [rPC_SM_train, varsd_train(:,sub_area{d})*v];
            rPC_SM_test = [rPC_SM_test, varsd_test(:,sub_area{d})*v];
        end            
    end

    rPC_BM_train = [];
    rPC_BM_test = [];
    for d=1:200
        temp = sum(Error_bm{d});
        N_bm(n,d) = find(temp==(min(temp)));

        BM_sub_train = NETd_node{d}(index_train,:);
        BM_sub_test = NETd_node{d}(index_test,:);

        netd_sub_train=nets_demean(BM_sub_train); % Gaussianise
        netd_sub_train = netd_sub_train/std(netd_sub_train(:));
        netd_sub_train = nets_demean(netd_sub_train-conf_train*(pinv(conf_train)*netd_sub_train));
        
        netd_sub_test=nets_demean(BM_sub_test); % Gaussianise
        netd_sub_test = netd_sub_test/std(netd_sub_test(:));
        netd_sub_test = nets_demean(netd_sub_test-conf_test*(pinv(conf_test)*netd_sub_test));
    
        COV = cov(netd_sub_train);
        [v,~] = eigs(COV, N_bm(n,d));
        if N_bm(n,d) ~= 1
%             v_r = rotatefactors(v);
%             rPC_BM_train = [rPC_BM_train, netd_sub_train*v_r];
%             rPC_BM_test = [rPC_BM_test, netd_sub_test*v_r];
            rPC_BM_train = [rPC_BM_train, netd_sub_train*v];
            rPC_BM_test = [rPC_BM_test, netd_sub_test*v];
        else
            rPC_BM_train = [rPC_BM_train, netd_sub_train*v];
            rPC_BM_test = [rPC_BM_test, netd_sub_test*v];
        end
    end

    [~, score] = pca(rPC_BM_train);
    rPC_BM_train100 = score(:,1:100);
    
    [~, score_test] = pca(rPC_BM_test);
    rPC_BM_test100 = score_test(:,1:100);
    %%% CCA
    [A_train{n}, B_train{n}, R_train{n}, P_train{n}, Q_train{n}] = canoncorr(rPC_SM_train, rPC_BM_train100);
    %%% reconstruct test canonical variates
    P_test{n} = rPC_SM_test*A_train{n};
    Q_test{n} = rPC_BM_test100*B_train{n};
    R_test{n} = corr(P_test{n}, Q_test{n});
    R_test{n} = diag(R_test{n});
    
    %%% permutation testing
    EB_in=EB(index_train,:);
    EB_out=EB(index_test,:);
    PAPset_in=palm_quickperms([ ], EB_in, Nperm);      % the final matrix of permuations
    PAPset_out=palm_quickperms([ ], EB_out, Nperm);     % the final matrix of permuations
    
    grotRp=[]; clear grotRpval grotRpval_out;
    grotRp_out=[];
    grotRp_cv=[];
    
    for j=1:Nperm
      [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(rPC_SM_train,rPC_BM_train100(PAPset_in(:,j),:));
      %Permutation for CV
       grotVr_CV=rPC_BM_test100(PAPset_out(:,j),:)*grotBr; %'psuedo-canonical variates' for SM
       grotUr_CV=rPC_SM_test*grotAr;
       for q=1:size(grotUr_CV,2)
            grotRp_cv(j,q)=corr(grotVr_CV(:,q), grotUr_CV(:,q));  %(1)
       end  
    end
    
    for i=1:size(grotRp,2);  % get FWE-corrected pvalues
      grotRpval(i)=(1+sum(grotRp(2:end,1)>=R_train{n}(i)))/Nperm;
      grotRpval_cv(i)=(1+sum(max(grotRp_cv(2:end,:)')>=R_test{n}(i)))/Nperm;  
    end
    Ncca_in(n)=sum(grotRpval<0.05)  % number of FWE-significant CCA components
    Ncca_cv(n)=sum(grotRpval_cv<0.05)  % number of FWE-significant CCA components
  
    %%% canonical loadings/structral coeff.
    CL_SM_train{n}=corr(P_train{n}(:,1:3),varsd_train,'rows','pairwise'); %in & in
    CL_BM_train{n}=corr(Q_train{n}(:,1:3),netd_train); %NETd(index1,1:size(NET,2))
    
    %%% Variance explained 
    TVE_SM_train=CL_SM_train{n}.^2; %in & in
    TVE_BM_train=CL_BM_train{n}.^2; %NETd(index1,1:size(NET,2))
   
    TVE_SM_test=[];
    TVE_BM_test=[]; 
    for j=1:3
        TVE_SM_test(j,:)=partialcorr(P_test{n}(:,j),varsd_test,P_test{n}(:,[1:j-1, j+1:3]),'rows','pairwise').^2; %in & in
        TVE_BM_test(j,:)=partialcorr(Q_test{n}(:,j),net_test,Q_test{n}(:,[1:j-1, j+1:3])).^2; %NETd(index1,1:size(NET,2))   
    end
    
    AVE_SM_train(n,:)=mean(TVE_SM_train,2,'omitnan');
    AVE_BM_train(n,:)=mean(TVE_SM_train,2,'omitnan');
     
    AVE_SM_test(n,:)=mean(TVE_SM_test,2,'omitnan');
    AVE_BM_test(n,:)=mean(TVE_SM_test,2,'omitnan');

end
