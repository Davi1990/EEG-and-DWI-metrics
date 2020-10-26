%%
%%Extract Network Values from matrices
%%DEFAULT%%
clear all
ROIs=importdata('Schaefer2018_ROIs_order.txt');
ROIs=ROIs.textdata
ROIs=ROIs(:,2);


for oo=1:22
    contenuto=importdata('lista.txt');
    matrix=char(contenuto(oo));
    matrix=importdata(matrix);
    %%
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
    %%
    %Connection rest of Brain Network
    Dafault_rest=Default_matrix(:,[1:37,51:89]);
    big=size(Dafault_rest);
    big=big(1,1)*big(1,2);
    s=sum(sum(Dafault_rest));
    CI_Default_between=s/big;
    %%
    %Connection with the other Networks
    Dafault_Vis=Default_matrix(:,[1:9,51:58]);
    Dafault_SomMot=Default_matrix(:,[10:15,59:66]);
    Dafault_DAN=Default_matrix(:,[16:23,67:73]);
    Dafault_SalVentAttn=Default_matrix(:,[24:30,74:78]);
    Dafault_Limbic=Default_matrix(:,[31:33,79:80]);
    Dafault_FPN=Default_matrix(:,[34:37,81:89]);
    %Dafault_Vis
    big=size(Dafault_Vis);
    big=big(1,1)*big(1,2);
    s=sum(sum(Dafault_Vis));
    CI_Dafault_Vis=s/big;
    %Dafault_SomMot
    big=size(Dafault_SomMot);
    big=big(1,1)*big(1,2);
    s=sum(sum(Dafault_SomMot));
    CI_Dafault_SomMot=s/big;
    %Dafault_DAN
    big=size(Dafault_DAN);
    big=big(1,1)*big(1,2);
    s=sum(sum(Dafault_DAN));
    CI_Dafault_DAN=s/big;
    %Dafault_SalVentAttn
    big=size(Dafault_SalVentAttn);
    big=big(1,1)*big(1,2);
    s=sum(sum(Dafault_SalVentAttn));
    CI_Dafault_SalVentAttn=s/big;
    %Dafault_Limbic
    big=size(Dafault_Limbic);
    big=big(1,1)*big(1,2);
    s=sum(sum(Dafault_Limbic));
    CI_Dafault_Limbic=s/big;
    %Dafault_FPN
    big=size(Dafault_FPN);
    big=big(1,1)*big(1,2);
    s=sum(sum(Dafault_FPN));
    CI_Dafault_FPNs=s/big;
    DMN_vs_other_6_Net = [CI_Dafault_Vis,CI_Dafault_SomMot,CI_Dafault_DAN,CI_Dafault_SalVentAttn,CI_Dafault_Limbic,CI_Dafault_FPNs]
    %%
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
    %%
    CI_ALL(oo,:)=[CI_Stimulation_complete, CI_Stimulation_DMN, CI_Stimulation_rest, CI_Default_within,CI_Default_between,CI_Dafault_Vis,CI_Dafault_SomMot,CI_Dafault_DAN,CI_Dafault_SalVentAttn,CI_Dafault_Limbic,CI_Dafault_FPNs]

end

%%

oo=1
contenuto=importdata('lista.txt');
matrix=char(contenuto(oo));
matrix=importdata(matrix);
myColorMap = lines(length(matrix));

ROIs=importdata('Schaefer2018_ROIs_order.txt');
ROIs=ROIs.textdata
ROIs=ROIs(:,2);
myLabel = ROIs
circularGraph(matrix,'Colormap',myColorMap,'Label',myLabel);



%%

%%Extract Network Values from matrices
%%DAN%%
clear all
ROIs=importdata('Schaefer2018_ROIs_order.txt');
ROIs=ROIs.textdata
ROIs=ROIs(:,2);


for oo=1:22
    contenuto=importdata('lista.txt');
    matrix=char(contenuto(oo));
    matrix=importdata(matrix);
    %%
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
    %%
    %Connection rest of Brain Network
    DAN_rest=DAN_matrix(:,[1:15,24:66,74:100]);
    big=size(DAN_rest);
    big=big(1,1)*big(1,2);
    s=sum(sum(DAN_rest));
    CI_DAN_between=s/big;
    %%
    %Connection with the other Networks
    DAN_Vis=DAN_matrix(:,[1:9,51:58]);
    DAN_SomMot=DAN_matrix(:,[10:15,59:66]);
    DAN_SalVentAttn=DAN_matrix(:,[24:30,74:78]);
    DAN_Limbic=DAN_matrix(:,[31:33,79:80]);
    DAN_FPN=DAN_matrix(:,[34:37,81:89]);
    DAN_Default=DAN_matrix(:,[38:50,90:100]);

    %DAN_Vis
    big=size(DAN_Vis);
    big=big(1,1)*big(1,2);
    s=sum(sum(DAN_Vis));
    CI_DAN_Vis=s/big;
    %DAN_SomMot
    big=size(DAN_SomMot);
    big=big(1,1)*big(1,2);
    s=sum(sum(DAN_SomMot));
    CI_DAN_SomMot=s/big;
    %DAN_SalVentAttn
    big=size(DAN_SalVentAttn);
    big=big(1,1)*big(1,2);
    s=sum(sum(DAN_SalVentAttn));
    CI_DAN_SalVentAttn=s/big;
    %DAN_Limbic
    big=size(DAN_Limbic);
    big=big(1,1)*big(1,2);
    s=sum(sum(DAN_Limbic));
    CI_DAN_Limbic=s/big;
    %DAN_FPN
    big=size(DAN_FPN);
    big=big(1,1)*big(1,2);
    s=sum(sum(DAN_FPN));
    CI_DAN_FPN=s/big;
    %DAN_Default
    big=size(DAN_Default);
    big=big(1,1)*big(1,2);
    s=sum(sum(DAN_Default));
    CI_DAN_Default=s/big;




    DAN_vs_other_6_Net = [CI_DAN_Vis,CI_DAN_SomMot,CI_DAN_SalVentAttn,CI_DAN_Limbic,CI_DAN_FPN,CI_DAN_Default]
    %%
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
    %%
    CI_ALL(oo,:)=[CI_Stimulation_complete, CI_Stimulation_DAN, CI_Stimulation_rest, CI_DAN_within,CI_DAN_between,CI_DAN_Vis,CI_DAN_SomMot,CI_DAN_SalVentAttn,CI_DAN_Limbic,CI_DAN_FPN,CI_DAN_Default]

end
