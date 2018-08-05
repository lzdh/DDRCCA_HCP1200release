clear 
clc
unrestricted = readtable('../rawdata/behaviour.csv');% download unrestricted behavioural data and save as 'behaviour.csv'
restricted  = readtable('../rawdata/restricted.txt');% download restricted behavioural data and save as 'behaviour.csv'

%%% merger restricted and unrestricted behavioural data
unrestricted.Subject = [];
behaviour_full = [restricted, unrestricted];

%%% Variable grouping
confounds = [202:204, 2,10,11,19:21, 223, 392,393];
Demographics = [12:18];
PhyHeal = [22,23,25:36];
FemHeal = [37:46];
FamHist = [47:64];
Psychiatry = [65:108];
Sensory = [109:113, 773:782];
Drug = [114:121,192:201];
Alcohol = [122:155];
Tobacco = [156:191];
Alertness = [288:315];
Cognition = [316:367];
Emotion = [368:391];
FS = [392:590];
Motor = [701:707];
Personality = [708:772];

%%%% construct and numeralise family info
Fam_info = behaviour_full(:,[1,4:8]);
Subtable=Fam_info.ZygosityGT;
A=zeros(1,size(Fam_info,1))'; % default not twin
index = find(strcmp(Subtable, 'DZ'));% i.e. DZ
A(index)=1;
index = find(strcmp(Subtable, 'MZ')); 
A(index)=2;
Fam_info.ZygosityGT=A;
%writetable(Fam_info,'Fam_info.csv');

%%%% select domains of interests
behaviour = behaviour_full(:,[1,confounds, Demographics,PhyHeal,FemHeal,FamHist,Psychiatry, ...
            Sensory,Drug,Alcohol,Tobacco,Alertness,Cognition,Emotion,Motor,Personality]);

%%%%%%%%%%conver into numerical matrix
 %%%Release
Subtable=behaviour.Release;
A=zeros(1,size(behaviour,1))';
  index = find(strcmp(Subtable, 'Q1'));
A(index)=1;
index = find(strcmp(Subtable, 'Q2'));
A(index)=2;
index = find(strcmp(Subtable, 'Q3'));
A(index)=3;
index = find(strcmp(Subtable, 'S500'));
A(index)=4;
index = find(strcmp(Subtable, 'MEG2'));
A(index)=5;
 index = find(strcmp(Subtable, 'S900'));
A(index)=6;
 index = find(strcmp(Subtable, 'S1200'));
A(index)=7;
 behaviour.Release=A;
 
 %%%Acquisition
Subtable=behaviour.Acquisition;
A=zeros(1,size(behaviour,1))';
  index = find(strcmp(Subtable, 'Q01'));
A(index)=1;
index = find(strcmp(Subtable, 'Q02'));
A(index)=2;
index = find(strcmp(Subtable, 'Q03'));
A(index)=3;
index = find(strcmp(Subtable, 'Q04'));
A(index)=4;
 index = find(strcmp(Subtable, 'Q05'));
A(index)=5;
 index = find(strcmp(Subtable, 'Q06'));
A(index)=6;
index = find(strcmp(Subtable, 'Q07'));
A(index)=7;
index = find(strcmp(Subtable, 'Q08'));
A(index)=8;
index = find(strcmp(Subtable, 'Q09'));
A(index)=9;
index = find(strcmp(Subtable, 'Q10'));
A(index)=10;
index = find(strcmp(Subtable, 'Q11'));
A(index)=11;
index = find(strcmp(Subtable, 'Q12'));
A(index)=12;
index = find(strcmp(Subtable, 'Q13'));
A(index)=13;
 behaviour.Acquisition=A;
 
%%%Eye correction: unknow coded as '999'. turn '999' to NaN
behaviour.Correction(behaviour.Correction == 999) = nan(1,length(behaviour.Correction(behaviour.Correction == 999)));

%%%Gender
  gender = behaviour.Gender;
 index = find(strcmp(gender, 'M'));
A=ones(1,size(behaviour,1))';
A(index)=0;
behaviour.Gender=A;

%%% Reconstruction variables
Subtable=behaviour.fMRI_3T_ReconVrs;
A=nan(1,size(behaviour,1))';
  index = find(strcmp(Subtable, 'r177'));
A(index)=1;
index = find(strcmp(Subtable, 'r177 r227'));
A(index)=2;
index = find(strcmp(Subtable, 'r227'));
A(index)=3;
behaviour.fMRI_3T_ReconVrs = A;

%%% Race
Subtable = behaviour.Race;
A=nan(1,size(behaviour,1))';
index_race = find(strcmp(Subtable, 'White'));
A(index_race)=1;
index_race = find(strcmp(Subtable, 'Am. Indian/Alaskan Nat.'));
A(index_race)=2;
index_race = find(strcmp(Subtable, 'Asian/Nat. Hawaiian/Othr Pacific Is.'));
A(index_race)=3;
index_race = find(strcmp(Subtable, 'Black or African Am.'));
A(index_race)=4;
index_race = find(strcmp(Subtable, 'More than one'));
A(index_race)=5;
behaviour.Race=A;

%%% Ethnicity
Subtable = behaviour.Ethnicity;
A=zeros(1,size(behaviour,1))';
index = find(strcmp(Subtable, 'Hispanic/Latino'));
A(index)=1;
behaviour.Ethnicity=A;

%%% Color vision
Subtable = behaviour.Color_Vision;
A=nan(1,size(behaviour,1))';
index = find(strcmp(Subtable, 'NORMAL'));
A(index)=1;
index = find(strcmp(Subtable, 'DEUTAN'));
A(index)=2;
index= find(strcmp(Subtable, 'PROTAN'));
A(index)=3;
index = find(strcmp(Subtable, 'TRITAN'));
A(index)=4;
behaviour.Color_Vision=A;

%%%Eye
Subtable = behaviour.Eye;
A=nan(1,size(behaviour,1))';
index = find(strcmp(Subtable, 'B'));
A(index)=0;
index = find(strcmp(Subtable, 'R'));
A(index)=1;
behaviour.Eye=A;

%%% Binary variables: true/false vars (binary)
B = behaviour(:,122:129);
sub_header = B.Properties.VariableNames;
A = zeros(size(behaviour,1),8);
B = table2cell(B);
for j= 1:8
index = find(strcmp(B(:,j), 'true'));
A(index,j)=1;
end
A=array2table(A);
A.Properties.VariableNames = sub_header;
behaviour= [behaviour(:,1:121) A behaviour(:,130:end)];

%%% Personality sub questions
A=[];
B = behaviour(:,326:385);
A = nan(size(behaviour,1),size(B,2));
sub_header = B.Properties.VariableNames;
B = table2cell(B);
for j=1:size(B,2)
    index = find(strcmp(B(:,j), 'SD'));
    A(index,j)=1;
    index = find(strcmp(B(:,j), 'D'));
    A(index,j)=2;
    index= find(strcmp(B(:,j), 'N'));
    A(index,j)=3;
    index = find(strcmp(B(:,j), 'A'));
    A(index,j)=4;
    index = find(strcmp(B(:,j), 'SA'));
    A(index,j)=5;
end
A=array2table(A);
A.Properties.VariableNames = sub_header;
behaviour= [behaviour(:,1:325) A];


%%% Create dummy variables for Race, Vision, etc.
Race_white = zeros(size(behaviour,1),1);
index = find(behaviour.Race == 1);
Race_white(index) = 1;

Race_black = zeros(size(behaviour,1),1);
index = find(behaviour.Race == 4);
Race_black(index) = 1;

Race_other = zeros(size(behaviour,1),1);
index = find(behaviour.Race == 2);
index = [index;find(behaviour.Race == 3)];
index = [index;find(behaviour.Race == 5)];
Race_other(index) = 1;

%%% Add them to behavioural_full and behaviour
Race_dummy = [Race_white, Race_black, Race_other];
Race_dummy = array2table(Race_dummy);
Race_dummy.Properties.VariableNames = {'Race_white', 'Race_black', 'Race_other'};
behaviour_full = [behaviour_full,Race_dummy];
behaviour = [behaviour, Race_dummy];

%%% Vision dummy
Vision_normal = zeros(size(behaviour,1),1);
index = find(behaviour.Color_Vision == 1);
Vision_normal(index) = 1;

Vision_deutan = zeros(size(behaviour,1),1);
index = find(behaviour.Color_Vision == 2);
Vision_deutan(index) = 1;

Vision_protan = zeros(size(behaviour,1),1);
index = find(behaviour.Color_Vision == 3);
Vision_protan(index) = 1;

Vision_tritan = zeros(size(behaviour,1),1);
index = find(behaviour.Color_Vision == 4);
Vision_tritan(index) = 1;

Vision_dummy = [Vision_normal, Vision_deutan, Vision_protan, Vision_tritan];
Vision_dummy = array2table(Vision_dummy);
Vision_dummy.Properties.VariableNames = {'Vision_normal', 'Vision_deutan', 'Vision_protan', 'Vision_tritan'};
behaviour_full = [behaviour_full,Vision_dummy];
behaviour = [behaviour, Vision_dummy];


%%% calculating sleeping hours
% PSQI_SleepHour = 24*datenum(behaviour.PSQI_GetUpTime)-24*datenum(behaviour.PSQI_BedTime);
% for i = 1:1206
%     if PSQI_SleepHour(i)<0
%         PSQI_SleepHour(i) = 24+PSQI_SleepHour(i);
%     end
% end

%%% Remove menstrual explain
behaviour.Menstrual_Explain = [];
behaviour.PSQI_BedTime = [];
behaviour.PSQI_GetUpTime = [];

%%% New grouping
confounds = [2:5,383:385,7:13];
Demographics = [14:20];
PhyHeal = [21:34];
FemHeal = [35:43];
FamHist = [44:61];
Psychiatry = [62:105];
Sensory = [386:389, 107:120];
Drug = [121:138];
Alcohol = [139:172];
Tobacco = [173:208];
Alertness = [209:234];
Cognition = [235:286];
Emotion = [287:310];
Motor = [311:317];
Personality = [318:382];

behaviour = behaviour(:,[1,confounds, Demographics,PhyHeal,FemHeal,FamHist,Psychiatry, ...
            Sensory,Drug,Alcohol,Tobacco,Alertness,Cognition,Emotion,Motor,Personality]);

headers = behaviour.Properties.VariableNames';

headers(216:222) = {'PSQI_SleepQuality1', 'PSQI_SleepLatency', 'PSQI_SleepQuality2', 'PSQI_SleepDuration', 'PSQI_SleepDisturbance'...
                    'PSQI_SleepMeds', 'PSQI_DayDysfunction'};
behaviour_num = table2array(behaviour);
%writetable(cell2table(headers), 'headers_SM387.txt', 'WriteVariableNames',false);

%% load brain measures
BM = importdata('../rawdata/netmats2.txt');
BM_ID = importdata('../rawdata/subjectIDs.txt');

%% Select subjects with both measures
[~,idx]= intersect(behaviour_num(:,1),BM_ID);
behaviour_num = behaviour_num(idx,:);
behaviour_full = behaviour_full(idx,:);
Fam_info = Fam_info(idx,:);
writetable(Fam_info, 'Fam_infoBMsubs.csv');
writetable(behaviour_full, 'SM1003x786.csv');

