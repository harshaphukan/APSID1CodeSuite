clc
clear
A=LoadData(); % Load Data for Layer 9
B=LoadData(); % Load Data for Layer10
%% Part 1 
%(Euler angles with 30degree correction for phi2)
Euler1=A(21,10:12); % Euler angle for parent 
EulerT=A(121,10:12);% Euler ange for twin 
Euler2=B(63,10:12); % Euler angle for neighbor
STens1=A(21,4:9); % Stress Tensor for parent
STensT=A(121,4:9); % Stress Tensor Twin
STens2=B(63,4:9); % Stress Tensor for neighbor

%% Part 5
% Variant Identification Step
% User-defined function returns array with slip/twin system and mprime info
[HT]=MPrimeFunHJP(Euler1,EulerT,STens1,STensT,0.2,1e-6);
% Loop reads thru the slip system/mprime table and reports the T1 system in
% twin that has mprime~1 w.r.t T1 system in parent grain
for i=1:length(H)
    if (HT(i,1)>=34 && HT(i,1)<=39) && (HT(i,3)>=34 && HT(i,3)<=39) && HT(i,end)>0.999
        ActiveVar=HT(i,1);
    end
end
 
MPrimeFun(Euler1,EulerT,Euler2,STens1,[0,0.25])


HexFun(Euler1,STens1)
%%
%This part calculate Euclidean distance between the parent grain and
%%neighbors in same layer and layers above and below it
Euclidean9=zeros(length(A),1);% Euclidean distance between parent and other grains (same layer)
for i=1:length(A)
Euclidean9(i)=sqrt((A(21,1)-A(i,1))^2+(A(21,2)-A(i,2))^2);
end
Euclidean9_10=zeros(length(B),1);
for i=1:length(B)
Euclidean9_10(i)=sqrt((A(21,1)-B(i,1))^2+(A(21,2)-B(i,2))^2);
end

D=[A,Euclidean9];
D=sortrows(D,13);
SL=D(2:end,:); % List of grains in same layer closest to parent
E=[B,Euclidean9_10];
LA=sortrows(E,13); % List of grains in layer above closest to parent
%% Part 4: This parts prunes the list of nearest neighbors
ProxSL=[];
for i=1:length(SL)
    if SL(i,end)<=100
        ProxSL(i)=SL(i,1);% Closest neighbors in same layer 
    end
end
ProxLA=[];
for i=1:length(LA)
    if LA(i,end)<=100
        ProxLA(i)=LA(i,1);% Closest neighbors in layer above
    end
end

% Generate Vectors linking COM of twin to COM of neighbors
TwinCOM=A(21,2:3); % Coordinates of Twin COM
SLTwin=UnitVec(A(ProxSL,2:3),TwinCOM);
SLTwin=[ProxSL',SLTwin]; % Appended Grain number to unit vector list
LATwin=UnitVec(B(ProxLA,2:3),TwinCOM);
LATwin=[ProxLA',LATwin];% Appended Grain number to unit vector list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine  Slip/Twin Systems that are within 70degrees w.r.t Vector connecting 
%COM of Twin and Neighbor in concern! Takes One Layer @ a time!
% Same Layer First!
NSDSL=zeros(57,3,length(SLTwin)); % Stacked array to store Slip/Twin  vectors for each neighbor insame layer 
SSList=zeros(1,57,length(SLTwin)); % Array to store Slip/Twin System # that meets criteria
F={};%zeros(300,7,length(SLTwin)); % Array to store mprime SS Info
K={};
%LogicVec=zeros(300,7,length(SLTwin));
for i=1:length(SLTwin)
[sortmvSL,slpsysSL]=HexEval(A(i,10:12),A(i,4:9));
SNSL=SlipSysFun(sortmvSL,slpsysSL);
for j=1:57
NSDSL(j,:,i)=SNSL(j).f3;
if (acosd(dot(SLTwin(i,2:4),NSDSL(j,:,i)))>=0) && (acosd(dot(SLTwin(i,2:4),NSDSL(j,:,i)))<=70)
SSList(1,j,i)=SNSL(j).f1;


end

end
[H]=MPrimeFunHJP(Euler1,A(i,10:12),STens1,A(i,4:9),0.2,1e-6);
MP{i}=H;
Z=SSList(:,:,i);
Z=Z(Z~=0);
K{i}=Z;
for k=1:length(MP{i})
    LogicVec=ismember(MP{i}(:,3),K{i});
    if (LogicVec(k)==1) && (MP{i}(k,end)>0.85) && (MP{i}(k,1)==ActiveVar)
F{i}=MP{i}(k,:);
    end
end
    

end


% F is the final pruned cell array showing the compatible slip system
% meeting the criteria imposed. The process can be repeated for layers
% above and below.



        