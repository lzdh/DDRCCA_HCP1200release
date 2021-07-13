clear
clc

SMvars = importdata('../Data/Flipped1003x318.mat');
headers = importdata('../Data/header318_flipped.txt'); 
behaviour_full = importdata('../Data/SM1003x789.mat');
NET = importdata('../rawdata/netmats2.txt');
SM_sdr = importdata('../Data/SDR_SM67.mat');
BM_sdr = importdata('../Data/SDR_BM6867.mat');

Nkeep = 67; % number of components to keep in the PCA of SM
Nperm = 10000; % number of permutations

%%% construct confounds
confounds = [SMvars(:,[2:6,8, 16:18, 316:318]), behaviour_full.FS_BrainSeg_Vol.^(1/3), behaviour_full.FS_IntraCranial_Vol.^(1/3)]; 
conf=palm_inormal(confounds); % Gaussianise.
conf(isnan(conf))=0;  % impute missing data as zeros
conf=nets_normalise([conf conf(:,7:end).^2]);  % add on squared terms and renormalise

%%% deconfound and normalise original SM and BM
varsd=palm_inormal(SMvars);
for i=1:size(varsd,2)
  grot=(isnan(varsd(:,i))==0); grotconf=nets_demean(conf(grot,:)); varsd(grot,i)=nets_normalise(varsd(grot,i)-grotconf*(pinv(grotconf)*varsd(grot,i)));
end

NET1=nets_normalise(NET); % no norm
NETd=nets_demean(NET1-conf*(pinv(conf)*NET1));   % deconfound and demean

%%% SVD/PCA of BM and SM
[uu1,ss1,~]=nets_svds(NETd,100); % SVD reduction
BM_pc = uu1*ss1;

COV = cov(varsd, 'partialrows'); % PCA
vv2 = eigs(COV, Nkeep);
SM_pc = varsd*vv2; % Principal components

%%% CCA
[A,B,R,P,Q] = canoncorr(SM_pc, BM_pc);

%%% CCA permutation testing
grotRp=zeros(Nperm,Nkeep); clear grotRpval;
for j=1:Nperm
  j
  [grotAr,grotBr,grotRp(j,:),grotUr,grotVr,grotstatsr]=canoncorr(SM_pc(PAPset(:,j),:), BM_pc);
end
for i=1:Nkeep;  % get FWE-corrected pvalues
  grotRpval(i)=(1+sum(grotRp(2:end,1)>=grotR(i)))/Nperm;
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


