clc
clear
A=LoadData(); % Load Data for Layer in which the twin/parent lies(1)
B=LoadData(); % Load Data for Layer below(2)
C=LoadData();

%%

  
for i=1:length(A)
[FA{i},~]=TwinMisor(A,A,A(i,1));

 [R1,C1]=find(FA{i}(:,2)>=83 & FA{i}(:,2)<=100);
  FA{i}=FA{i}(R1,:);

  [~,FB{i}]=TwinMisor(A,B,A(i,1));

 [R2,C2]=find(FB{i}(:,2)>=83 & FB{i}(:,2)<=100);
  FB{i}=FB{i}(R2,:);
  
  [~,FC{i}]=TwinMisor(A,C,A(i,1));

 [R3,C3]=find(FC{i}(:,2)>=83 & FC{i}(:,2)<=100);
  FC{i}=FC{i}(R3,:);
  

end


%%
GrainPlotter(A(77,2:3),C(28,2:3),77,28)
%%
MPrimeFun([0,0,0],A(77,10:12),[0,0,0],C(1,10:12),C(1,4:9),[0,0.25],6)

%%
HexFunRot(A(85,10:12),A(85,4:9))
%HexFun(A(15,10:12),A(15,4:9))
