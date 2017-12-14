% Code to systematically analyze compatible slip systems between identified
% parent grains and neighbors
% Has capability of analyzing two layers at a time
clc
clear
A=LoadData(); % L(3)
%B=LoadData(); % Load Data for Layer below(2)
%C=LoadData(); %L(4)
%%
%[vMArray]=vonMises(A)

%% Part 1 
%(Euler angles with 30degree correction for phi2)
PIndex=3; % Grain number of Parent 
TIndex=32; % Grain number of Twin
Euler1=A(PIndex,10:12); % Euler angle for parent 
EulerT=A(TIndex,10:12);% Euler ange for twin 
STens1=A(PIndex,4:9); % Stress Tensor for parent
STensT=A(TIndex,4:9); % Stress Tensor Twin
%%
[OLA,OLB]=TwinMisor(A,A,TIndex)
%% Part 5
% Variant Identification Step
% User-defined function returns array with slip/twin system and mprime info
[HT]=MPrimeFunHJP(Euler1,EulerT,STens1,STensT,0.2,1e-6);
% Loop reads thru the slip system/mprime table and reports the T1 system in
% twin that has mprime~1 w.r.t T1 system in parent grain
for i=1:length(HT)
    if (HT(i,1)>=34 && HT(i,1)<=39) && (HT(i,3)>=34 && HT(i,3)<=39) && HT(i,end)>0.999
        ActiveVar=HT(i,1);
    end
end
 

%%
MPrimeFun(Euler1,EulerT,[0,0,0],[0,0,0],STens1,[0,0.25],2)
%%
HexFun(EulerT,STensT)
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
EuclideanC(i)=sqrt((B(PIndex,2)-C(i,2))^2+(B(PIndex,3)-C(i,3))^2);
end
%}
%{
D=[C,EuclideanC];
LB=sortrows(D,13);% List of grains in layer below(1)
%}

E=[B,EuclideanB];
LA=sortrows(E,13); % List of grains in layer above(2)

W=[A,EuclideanA];
W=sortrows(W,13); % List of grains in same layer as parent/twin(1)
SL=W(2:end,:);

TwinCOM=A(TIndex,2:3); % Coordinates of Twin COM

%% Part 4: This parts prunes the list of nearest neighbors
for i=1:length(SL)
    if SL(i,end)<=100
        ProxSL(i)=SL(i,1);% Closest neighbors in same layer(4)
    else
        ProxSL(i)=0;
    end
end
ProxSL=ProxSL(ProxSL~=0);

for i=1:length(LA)
    if LA(i,end)<=100
        ProxLA(i)=LA(i,1);% Closest neighbors in layer above (5)
    else
        ProxLA(i)=0;
    end
end
ProxLA=ProxLA(ProxLA~=0);

%{
for i=1:length(LB)
    if LB(i,end)<=100
        ProxLB(i)=LB(i,1);% Closest neighbors in layer below (3)
    end
end
ProxLB=ProxLB(ProxLB~=0);
%}
%%
% Generate Vectors linking COM of twin to COM of neighbors
if isempty(ProxSL)~=1
LogicV=find(ProxSL~=TIndex); % Ensure that the twin is not included here
ProxSL=ProxSL(LogicV);
SLTwin=UnitVec(A(ProxSL,2:3),TwinCOM);
SLTwin=[ProxSL',SLTwin]; % Appended Grain number to unit vector list
end

if isempty(ProxLA)~=1
LATwin=UnitVec(B(ProxLA,2:3),TwinCOM);
LATwin=[ProxLA',LATwin];% Appended Grain number to unit vector list
end
%{
if isempty(ProxLB)~=1
LBTwin=UnitVec(B(ProxLB,2:3),TwinCOM);
LBTwin=[ProxLB',LBTwin];% Appended Grain number to unit vector list
end
%}
%%
NIndexList=LATwin(:,1);% Grain #s of nearest neighbors shortlist
for i=1:length(NIndexList)
    FF{i}=SSPruner(B,Euler1,STens1,ActiveVar,NIndexList(i));
end

%%
% Function to tabulate slip system data for each shortlisted neighbor (in
% this for the same layer)
ii=6;
TabulateMPrime(FF{ii})
%%
GrainPlotter(A(49,2:3),A(124,2:3),49,124)
%% 1,2,3,4

MPrimeFun(Euler1,[0,0,0],[0,0,0],A(NIndexList(ii,1),10:12),STens1,[0,0.25],4)
%%
HexFunRot(A(24,10:12),A(24,4:9))