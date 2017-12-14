%Creates data structure containing slip system no,plane normal, slip
%direction & Schmid Factor using output from HexEval.m
function[S]=SlipSysFun(sortmv,slpsys) 
nslphex=57;
b_hkil=zeros(nslphex,4);
b=zeros(nslphex,3);
n_hkil=zeros(nslphex,4);
n=zeros(nslphex,3);
c_a=1.587;
% Convert Four Index Notation to Cartesian Coordinates!
for u=1:1:nslphex
    n_hkil(u,1)=slpsys{2,u}(1,1); n_hkil(u,2)=slpsys{2,u}(1,2);
    n_hkil(u,4)=slpsys{2,u}(1,4);
    n(u,1)=n_hkil(u,1);%2*n_hkil(u,1)+n_hkil(u,2)
    n(u,2)=(n_hkil(u,1)+2*n_hkil(u,2))/sqrt(3);%n_hkil(u,1)+2*n_hkil(u,2)
    n(u,3)=n_hkil(u,4)/c_a;%n_hkil(u,4)
    %% Change from here 11/29/2016
    b_hkil(u,1)=slpsys{2,u}(2,1); b_hkil(u,2)=slpsys{2,u}(2,2);
    b_hkil(u,4)=slpsys{2,u}(2,4);

    b(u,1)= b_hkil(u,1);%2*b_hkil(u,1)+b_hkil(u,2)  
    b(u,2)=(b_hkil(u,1)+2*b_hkil(u,2))/sqrt(3);%b_hkil(u,1)+2*b_hkil(u,2)   
    b(u,3)=b_hkil(u,4)/c_a;%b_hkil(u,4)
    %% End change 11/29/2016
    b(u,:)=slpsys{3,u}(2,:);
end

sort_Schm=sortmv(1:57,:);
Schm=sortrows(sort_Schm,1);

f1='Slip_system_no';
f2='Plane_Normal';f3='Burger_Vector';f4='Schmid_Factor';
v1=zeros(nslphex,1); v2=zeros(nslphex,3); v3=zeros(nslphex,3);v4=zeros(nslphex,1);

S=struct(f1,v1,f2,v2,f3,v3,f4,v4); %Data Structure to store Slip System information and Schmid Factor for a Grain
for u=1:1:nslphex
    v1=u; v2=n(u,:)/norm(n(u,:)); %Convert Plane normal into unit vector
    v3=b(u,:);%/norm(b(u,:));%Convert Slip direction/Burger's vector into unit vector
    v4=Schm(u,3);
    S(u).f1=v1;  %Slip system number
    S(u).f2=v2;  %Slip plane normal
    S(u).f3=v3;  %Slip direction
    S(u).f4=v4;  %Schmid Factor
end

end
