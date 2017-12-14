%Function to identify Twin Parent Pairs
%Inputs: Arrays from LoadData Function
%A-->Layer where twin is identified
%B--> Layer above (or below) or same layer as A
%Col1: Grain #in Layer, Col2,3:X,Y Coords, Cols 4:9:6 component of Stress
%Tensor, Col 10:12: Bunge Euler angles

function [OLA,OLB]=TwinMisor(A,B,TIndex)
%This part calculate Euclidean distance between the parent grain and
%%neighbors in same layer and layers above and below it
EuclideanA=zeros(length(A),1);% Euclidean distance between parent and other grains (Layer )
for i=1:length(A)
EuclideanA(i)=sqrt((A(TIndex,2)-A(i,2))^2+(A(TIndex,3)-A(i,3))^2);
end

EuclideanB=zeros(length(B),1);
for i=1:length(B)
EuclideanB(i)=sqrt((A(TIndex,2)-B(i,2))^2+(A(TIndex,3)-B(i,3))^2);
end

E=[B,EuclideanB];
LB=sortrows(E,13); 

W=[A,EuclideanA];
W=sortrows(W,13); 
LA=W(2:end,:);

LVA=find(LA(LA(:,end)<=250));
LVB=find(LB(LB(:,end)<=250));
LA=LA(LVA,:);
LB=LB(LVB,:);

OLA=zeros(length(LVA),6);
for i=1:length(LVA)
    OLA(i,1)=LA(i,1);
    [OLA(i,2),OLA(i,3),OLA(i,4),OLA(i,5),OLA(i,6)]=c_a_CalcFun(A(TIndex,10:12),LA(i,10:12));
end
OLB=zeros(length(LVB),6);
for i=1:length(LVB)
    OLB(i,1)=LB(i,1);
    [OLB(i,2),OLB(i,3),OLB(i,4),OLB(i,5),OLB(i,6)]=c_a_CalcFun(A(TIndex,10:12),LB(i,10:12));
end





end