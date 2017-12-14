%clear
clc

load('Layer02HMPLayer04GrainNum.mat') 

load('ID1StrainInfo.mat');
load('LogFileListL4.mat');

AllStrain=ID1StrainInfo(:,2);

Strain=[0.,0.03,0.03,0.07,0.14,0.22,1.28,1.48,1.56,2.04,2.17,2.42,2.50,2.60,2.93];

Index=find(ismember(AllStrain,Strain));

GNum=GrainNum(Index)

LF=LogFile(Index)

%{
    'CPTi3_0099_t100_cut_7.log'
    'CPTi3_00110_t100_cut_7.log'
    'CPTi3_00121_t100_cut_7.log'
    'CPTi3_00136_t100_cut_7.log'
    'CPTi3_00147_t100_cut_7.log'
    'CPTi3_00335_t100_cut_6.log'
    'CPTi3_00390_t100_cut_6.log'
    'CPTi3_00401_t100_cut_5.log'
    'CPTi3_00412_t100_cut_5.log'
    'CPTi3_00425_t100_cut_5.log'
    'CPTi3_00436_t100_cut_5.log'
    
%}

