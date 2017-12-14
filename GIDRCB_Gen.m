% Function to generate RCB and Grain ID arrays given the Euler angles as
% inputs
function [IDgr,RBdy]=GIDRCB_Gen(EulerTL,EulerTR,EulerBL,EulerBR)
%       ^
%  1    |  2
%<-------------->
%  3    |  4


EulerTLR=deg2rad(EulerTL);
EulerTRR=(pi/180)*EulerTR;
EulerBLR=(pi/180)*EulerBL;
EulerBRR=(pi/180)*EulerBR;


IDgr=zeros(4,11);

IDgr(:,1)=1:4;
IDgr(1,2:4)=EulerTL;
IDgr(2,2:4)=EulerTR;
IDgr(3,2:4)=EulerBL;
IDgr(4,2:4)=EulerBR;
IDgr(:,5)=[10 20 10 20];
IDgr(:,6)=[10 10 20 20];
IDgr(:,7:9)=rand(4,3);
IDgr(:,11)=[1 1 0 0];

RBdy=ones(6,14);
RBdy(1,1:3)=EulerTRR;
RBdy(1,4:6)=EulerTLR;
RBdy(2,1:3)=EulerBRR;
RBdy(2,4:6)=EulerBLR;
RBdy(3,1:3)=EulerBLR;
RBdy(3,4:6)=rand(1,3);
RBdy(4,4:6)=RBdy(3,4:6);
RBdy(4,1:3)=rand(1,3);
RBdy(5,1:3)=RBdy(4,1:3);
RBdy(6,1:3)=RBdy(2,4:6);
RBdy(6,4:6)=rand(1,3);
RBdy(:,7)=2*RBdy(:,7);
RBdy(:,8)=[0,90,135,45,90,0];
RBdy(:,9)=[5,15,10,10,15,15];
RBdy(:,10)=[15,5,10,20,15,15];
RBdy(:,11)=[15,15,20,20,15,25];
RBdy(:,12)=[15,15,20,10,25,15];
RBdy(:,13)=[3,2,2,4,4,4];
RBdy(:,14)=[1,1,3,1,3,2];
end