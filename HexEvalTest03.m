% Code to systematically analyze compatible slip systems between identified
% parent grains and neighbors
% Has capability of analyzing two layers at a time
clc
%clear
FPA='/Users/harshaphukan/Documents/Research/BeamLine1Scripts/Layer1/';
LFA='CPTi3_00315_t100_cut_7.log';

%{
FPB='/Users/harshaphukan/Documents/Research/BeamLine1Scripts/Layer11/';
LFB='CPTi3_00100_t100_cut_7.log';

%FPC='/Users/harshaphukan/Documents/Research/BeamLine1Scripts/Layer05/';
%LFC='CPTi3_00319_t100_cut_7.log';
%}
%B=LoadDataSameDir(FPB,LFB); % Load Data for Layer Above(2)
A=LoadDataSameDir(FPA,LFA);
%C=LoadDataSameDir(FPC,LFC);

%% Part 1 
%(Euler angles with 30degree correction for phi2)
PIndex=12; % Grain number of Parent 
TIndex=121; % Grain number of Twin
Euler1=A(PIndex,10:12); % Euler angle for parent 
EulerT=A(TIndex,10:12);% Euler ange for twin 
STens1=A(PIndex,4:9); % Stress Tensor for parent
STensT=A(TIndex,4:9); % Stress Tensor Twin
%% Part 5
% Variant Identification Step
% User-defined function returns array with slip/twin system and mprime info
[HT]=MPrimeFunHJP(Euler1,EulerT,STens1,STensT,0.2,1e-6);
% Loop reads thru the slip system/mprime table and reports the T1 system in
% twin that has mprime~1 w.r.t T1 system in parent grain
for i=1:length(HT)
    if (HT(i,1)>=34 && HT(i,1)<=39) && (HT(i,3)>=34 && HT(i,3)<=39) && HT(i,end)>=0.999
        ActiveVar=HT(i,1);
    end
end
 

%%
MPrimeFun(Euler1,EulerT,[0,0,0],[0,0,0],STens1,[0,0.25],2)
%%
HexFunRot(EulerT,STensT)
%%
%This part calculate Euclidean distance between the parent grain and
%%neighbors in same layer and layers above and below it
EuclideanA=zeros(length(A),1);% Euclidean distance between parent and other grains (Layer )
for i=1:length(A)
EuclideanA(i)=sqrt((A(PIndex,2)-A(i,2))^2+(A(PIndex,3)-A(i,3))^2);
end

EuclideanB=zeros(length(B),1);
for i=1:length(B)
EuclideanB(i)=sqrt((A(PIndex,2)-B(i,2))^2+(A(PIndex,3)-B(i,3))^2);
end



%{
EuclideanC=zeros(length(C),1);
for i=1:length(C)
EuclideanC(i)=sqrt((A(PIndex,2)-C(i,2))^2+(A(PIndex,3)-C(i,3))^2);
end


D=[C,EuclideanC];
D=sortrows(D,13);
LA=D(1:end,:);% List of grains in layer above

%}

E=[B,EuclideanB];
E=sortrows(E,13);
LB=E(1:end,:); % List of grains in layer above(2)

W=[A,EuclideanA];
W=sortrows(W,13); % List of grains in same layer as parent/twin(1)
SL=W(2:end,:);

TwinCOM=A(TIndex,2:3); % Coordinates of Twin COM


%% Simpler version of pruning the list of nearest neigbors in the 3 layers

[IVSL,~]=find(SL(:,13)<=120);
[IVLB,~]=find(LB(:,13)<=120);
%[IVLA,~]=find(LA(:,13)<=120);

ProxSL=SL(IVSL,1);
%ProxLA=LA(IVLA,1);
ProxLB=LB(IVLB,1);


%%
% Generate Vectors linking COM of twin to COM of neighbors
if isempty(ProxSL)~=1
LogicV=find(ProxSL~=TIndex); % Ensure that the twin is not included here
ProxSL=ProxSL(LogicV);
SLTwin=UnitVec(A(ProxSL,2:3),TwinCOM);
SLTwin=[ProxSL,SLTwin]; % Appended Grain number to unit vector list
end

if isempty(ProxLA)~=1
LATwin=UnitVec(B(ProxLA,2:3),TwinCOM);
LATwin=[ProxLA,LATwin];% Appended Grain number to unit vector list
end

if isempty(ProxLB)~=1
LBTwin=UnitVec(B(ProxLB,2:3),TwinCOM);
LBTwin=[ProxLB,LBTwin];% Appended Grain number to unit vector list
end

%%
NIndexList=ProxLB;% Grain #s of nearest neighbors shortlist
for i=1:length(NIndexList)
    FF{i}=SSPruner(B,Euler1,STens1,ActiveVar,NIndexList(i));
end

%%
% Function to tabulate slip system data for each shortlisted neighbor (in
% this for the same layer)
ii=4;
TabulateMPrime(FF{ii})
%%
GrainPlotter(B(17,2:3),A(1,2:3),17,1)
%%
GrainPlotter(A3(57,2:3),A4(7,2:3),57,7)
%% 1,2,3,4
Z1=[0,0,0];

MPrimeFun(Euler1,Z1,Z1,B(NIndexList(ii,1),10:12),STens1,[0,0.25],4)
%% 
%%
HexFunRot(C(NIndexList(ii,1),10:12),C(NIndexList(ii,1),4:9))
%%
HexFunRot(B(12,10:12),B(12,4:9))