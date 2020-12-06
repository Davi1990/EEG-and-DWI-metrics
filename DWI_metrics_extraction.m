%%

%%Author: Davide Momi, PhD [momi.davide89@gmail.com],
%%https://twitter.com/davemomi
%%https://davi1990.github.io/

%for updatated version please refer to 
%https://github.com/Davi1990/DissNet

%%Extract Network Values from matrices
%%DEFAULT%%
clear all
addpath('BCT')
load('ROIs.mat');


load('data_exp.mat');
matrix = data_exp;
%
%Default
lh_Default=ROIs(38:50,:);
rh_Default=ROIs(90:100,:);
Default=[lh_Default;rh_Default];
lh_Default_matrix=matrix(38:50,:);
rh_Default_matrix=matrix(90:100,:);
Default_matrix=[lh_Default_matrix;rh_Default_matrix];
%Connection within DMN
Dafault_puro= Default_matrix(:,[38:50,90:100]);
big=size(Dafault_puro);
big=big(1,1);
big=(big*big)-big(1,1);
s=sum(sum(Dafault_puro));
CI_Default_within=s/big;
%

%Estract connectivity from stimulation point (rh_DMN_basically)
Stimulation_DMN= Dafault_puro(14,:);
big=size(Stimulation_DMN);
big=big(1,1)*big(1,2);
s=sum(sum(Stimulation_DMN));
CI_Stimulation_DMN=s/big;
Stimulation_rest=Default_matrix(14,[1:37,51:89]);
big=size(Stimulation_rest);
big=big(1,1)*big(1,2);
s=sum(sum(Stimulation_rest));
CI_Stimulation_rest=s/big;
CI_Stimulation_complete = matrix(90,:);
CI_Stimulation_complete = sum(CI_Stimulation_complete);
CI_Stimulation_complete = CI_Stimulation_complete/99
%

%Modularity
[Ci,Q]=modularity_und(matrix);
%
whole_brain = mean(matrix(:));

colnames = { 'Stim2Brain', 'Stim2Network', ...
            'Network', 'Modularity', 'Brain' };

        
 DMN_ALL= array2table([CI_Stimulation_complete, ...
                    CI_Stimulation_DMN, CI_Default_within, ...
                    Q, whole_brain], ...
                'VariableNames', colnames)



%


%%

%%Extract Network Values from matrices
%%DAN%%

clear all
load('ROIs.mat');


load('data_exp.mat');
matrix = data_exp;


%
%DAN
lh_DAN=ROIs(16:23,:);
rh_DAN=ROIs(67:73,:);
DAN=[lh_DAN;rh_DAN];
lh_DAN_matrix=matrix(16:23,:);
rh_DAN_matrix=matrix(67:73,:);
DAN_matrix=[lh_DAN_matrix;rh_DAN_matrix];
%Connection within DAN
DAN_puro= DAN_matrix(:,[16:23,67:73]);
big=size(DAN_puro);
big=big(1,1);
big=(big*big)-big(1,1);
s=sum(sum(DAN_puro));
CI_DAN_within=s/big;
%
%Connection rest of Brain Network
DAN_rest=DAN_matrix(:,[1:15,24:66,74:100]);
big=size(DAN_rest);
big=big(1,1)*big(1,2);
s=sum(sum(DAN_rest));
CI_DAN_between=s/big;
%
%Estract connectivity from stimulation point (RH_DorsAttn_Post_3)
Stimulation_DAN= DAN_puro(11,:);
big=size(Stimulation_DAN);
big=big(1,1)*big(1,2);
s=sum(sum(Stimulation_DAN));
CI_Stimulation_DAN=s/big;
Stimulation_rest=DAN_matrix(11,[1:15,24:66,74:100]);
big=size(Stimulation_rest);
big=big(1,1)*big(1,2);
s=sum(sum(Stimulation_rest));
CI_Stimulation_rest=s/big;
CI_Stimulation_complete = matrix(69,:);
CI_Stimulation_complete = sum(CI_Stimulation_complete);
CI_Stimulation_complete = CI_Stimulation_complete/99
%
%Modularity
[Ci,Q]=modularity_und(matrix);
%
whole_brain = mean(matrix(:));

colnames = { 'Stim2Brain', 'Stim2Network', ...
            'Network', 'Modularity', 'Brain' };

 DAN_ALL= array2table([CI_Stimulation_complete, ...
                    CI_Stimulation_DAN, CI_DAN_within, ...
                    Q, whole_brain], ...
                'VariableNames', colnames)

            
 