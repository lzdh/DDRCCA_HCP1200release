clear
clc
load('../Data/NETd_node.mat')
Fam_info=readtable('../Data/Fam_infoBMsubs.csv');
[C0,IA]=unique(Fam_info.Family_ID);
Subject_ID = importdata('../rawdata/subjectIDs.txt');

K=5; % K-fold CV
n=getenv('SGE_TASK_ID');

Node=NETd_node{n};

%%% normalise each node
Noded=nets_normalise(Node); 
Noded(isnan(Noded)) = 0;

%%% Calculate dim. (minimising  PRESS) for each of the sub-areas in held-in set using LOOCV!

for k=1:10

    Error2_sim = [];

    groups = crossvalind('Kfold',length(C0),K); %K-fold
    Fam_info.group_ID=zeros(1003,1);
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
            index2(i)=find(Test(i)==Subject_ID(:,1));
        end
        Noded_out=Noded(index2',:);

        for i=1:size(Train,1)
            index1(i)=find(Train(i)==Subject_ID(:,1));
        end
        Noded_in=Noded(index1',:);

        [~,dd1,vv1]=svd(Noded_in);
        if dd1(1,1)<dd1(2,2)  vv1=fliplr(vv1); end


        err2=[];
        for j=1:100             
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
Error2 = sum(Error2_ave); 
N_pca=find(Error2==(min(Error2)));
    


savedir='Outputs';
save(fullfile(savedir,['SDRNET_node' n]),'Error2','N_pca');

  
