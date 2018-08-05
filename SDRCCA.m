clear
clc

SMvars = importdata('../Data/Flipped1003x233.mat');
headers = importdata('../Data/header233_flipped.txt'); 
NET = importdata('../rawdata/netmats2.txt');
N_SM = importdata('../Data/N_pca_SM.mat');
N_BM = importdata('../Data/N_pca_BM.mat');
load('../Data/NETd_node.mat');
%SM_sdr = importdata('../Data/SDR_SM62.mat');
%BM_sdr = importdata('../Data/SDR_BM6867.mat');

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

%%% setup confounds matrix
conf=palm_inormal([SMvars(:,confounds(1:12)), SMvars(:,confounds(13:14)).^(1/3)]); 
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf, conf(:,5).^2, conf(:,10:12).^2]);  % add on squared terms and renormalise

sub_area={Demographics,PhyHeal,FemHeal, FamHist, Psychiatry, Sensory, Drug,Alcohol, Tobacco, Alertness, Cognition, Emotion, ...
             Motor, Personality};

rPC_SM = [];
rPC_SM_p = [];
for d=1:length(sub_area)
    SM_sub = SMvars(:,sub_area{d});
        
    varsd=palm_inormal(SM_sub); % Gaussianise
    for i=1:size(varsd,2) % deconfound ignoring missing data
      grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
    end
    
    [v_p,score,~]=ppca(varsd, N_SM(d));
    
    %%%impute missing data by k nearest rows
    if d == 3 %feminine health
        varsd(isnan(varsd))=0;
%         for j = 1:size(varsd,2)
%             varsd(isnan(varsd(:,j)),j)=mean(varsd(:,j),'omitnan');
%         end
    end
    varsd = knnimpute(varsd, 3);
    
    COV = cov(varsd, 'partialrows');
    [v,~] = eigs(COV, N_SM(d));
    
    varsd(isnan(varsd))=0;
    
    if N_SM ~= 1
        [v_r,T] = rotatefactors(v);
        rPC_SM = [rPC_SM, varsd*v_r];
        rPC_SM_p = [rPC_SM_p, score*T];
    else
        rPC_SM = [rPC_SM, varsd*v];
        rPC_SM_p = [rPC_SM_p, score];
    end
end

rPC_BM = [];
for d=1:200
    BM_sub = NETd_node{d};
        
    varsd=palm_inormal(BM_sub); % Gaussianise
    for i=1:size(varsd,2) % deconfound ignoring missing data
      grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
    end
    
    COV = cov(varsd, 'partialrows');
    [v,~] = eigs(COV, N_BM(d));
    if N_SM ~= 1
        v_r = rotatefactors(v);
        rPC_BM = [rPC_BM, varsd*v_r];
    else
        rPC_BM = [rPC_BM, varsd*v];
    end
end
    
Nkeep = 62; % number of components to keep in the PCA of SM
Nperm = 10000; % number of permutations

%%% deconfound and normalise original SM and BM
varsd=palm_inormal(SMvars);
for i=1:size(varsd,2)
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=nets_normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end

NET1=nets_demean(NET);  NET1=NET1/std(NET1(:)); % no norm
NETd=nets_demean(NET1-conf*(pinv(conf)*NET1));   % deconfound and demean

%%% De-confound SDR SM and BM
SM_sdr = rPC_SM-conf*(pinv(conf)*rPC_SM);
BM_sdr = rPC_BM-conf*(pinv(conf)*rPC_BM);

[~,BM_score] = pca(BM_sdr);
BM_pc = BM_score(:,1:100);

%%% CCA
[A,B,R,P,Q] = canoncorr(SM_sdr, BM_pc);

%%% CCA permutation testing

EB=hcp2blocks('../rawdata/restricted.txt', [ ], false, SMvars(:,1)); % change the filename to your version of the restricted file
PAPset=palm_quickperms([ ], EB, Nperm);                                            % the final matrix of permuations

grotRp=zeros(Nperm,Nkeep); clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(SM_sdr(PAPset(:,j),:), BM_pc);
end
for i=1:Nkeep;  % get FWE-corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=R(i)))/Nperm;
end
grotRpval
Ncca=sum(grotRpval<0.05)  % number of FWE-significant CCA components

%%% CCA loadings/structural coefficient
for i=1:3
    CL_bm(:,i) = corr(Q(:,i),NETd)';
    CL_sm(:,i) = corr(P(:,i),varsd,'rows','pairwise')'; % weights after deconfounding
end

%%% factor rotation
[CL_sm_r, T_sm] = rotatefactors(CL_sm);
[CL_bm_r, T_bm] = rotatefactors(CL_bm);



