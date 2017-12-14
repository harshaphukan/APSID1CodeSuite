function [HP]=MPrimeFunHJP(Euler1,Euler2,STens1,STens2,AvgSFTol,SFTol)
%Orientation Matrix
g1=OrMat(Euler1);
g2=OrMat(Euler2);

[sortmv1,slpsys1]=HexEval(Euler1,STens1);%This function obtains the Schmid matrices, stress tensors and euler angles using .log and .gve files as input

[S_A]=SlipSysFun(sortmv1,slpsys1);

[sortmv2,slpsys2]=HexEval(Euler2,STens2);%This function obtains the Schmid matrices, stress tensors and euler angles using .log and .gve files as input
[S_B]=SlipSysFun(sortmv2,slpsys2);

nslphex=57;
m_prime=cell(1,nslphex);%Declaring Array containing Slip System number, m-prime and average Schmid Factor
%Next Step is to evaluate m-prime values and tabulate according to average Schmid Factor
 for i=1:nslphex
     for j=1:nslphex
        m_prime{:,j}(i,1)=S_A(j).f1; %Slip System Number of Grain A
        m_prime{:,j}(i,2)=S_A(j).f4; %Schmid Factor for jth Slip System of Grain A
        m_prime{:,j}(i,3)=S_B(i).f1;  % ith Slip System Number of Grain B
        m_prime{:,j}(i,4)=S_B(i).f4;  %Schmid Factor for ith Slip System of Grain B
        m_prime{:,j}(i,5)=(S_A(j).f4 + S_B(i).f4)*0.5; %Calculation of Average Schmid Factor between jth slip system of A and ith slip system of B
        %m_prime{:,j}(i,6)=dot(log_data(Grain_indexA).U*S_A(j).f2',log_data(Grain_indexB).U*S_B(i).f2')*dot(log_data(Grain_indexA).U*S_A(j).f3',log_data(Grain_indexB).U*S_B(i).f3');
        %%m-prime calculation:jth ss of A and ith ss of B:Used Luster-Morris Convention
        m_prime{:,j}(i,6)=dot(S_A(j).f2*g1,S_B(i).f2*g2)*dot(S_A(j).f3*g1,S_B(i).f3*g2);
       
        
    end
end
%Next,we sort the m_prime matrices in decreasing order of average Schmid
%Factor
sort_mp=cell(1,nslphex); %Preallocating cell array for storing m-prime/slip system information in descending order of average Schmid Factor
for i=1:nslphex
    sort_mp{:,i}=sortrows(m_prime{:,i},-5);
end
%Input tolerence value for Average Schmid Factor

K=cell(1,nslphex); %Pre-allocate cell array to store slip system information based on a tolerance value of Average Schmid Factor




for i=1:nslphex
    for j=1:nslphex
        if sort_mp{:,j}(i,5)>=AvgSFTol && sort_mp{:,j}(i,6)~=0 %Checks Avg Schmid Factor with Schm_tol and ensures only rows with nonzero m-prime values are considered
           K{:,j}(i,:)=sort_mp{:,j}(i,:); %Populate slip system information for each of the 57 slip systems corresponding to the parent grain!
           K{:,j}(all(~K{:,j},2),:)=[]; %Remove all the zero rows
        end
    end
end

% Populate the Slip System data on the basis of Average Schmid Factor onto
% single matrix
H=K{:,1}; %Assign the array corresponding to first slip system of parent grain to the matrix H

ctr=2; %Counter to maintain track of the number of slip systems in the parent grain
while ctr<=57   %Total number of slip/twin systems in the parent grain
   
      H=[H;K{:,ctr}]; %Append H
      ctr=ctr+1;
end

L=size(H);

%Basal and First order pyramidal slip system information here
H1stOrder=zeros(100,6);
for i=1:L(1)
    if (H(i,1)==1 || H(i,1)==2 || H(i,1)==3 || H(i,1)==4 ||  H(i,1)==5 || H(i,1)==6) && (H(i,3)==1 || H(i,3)==2 ...
            || H(i,3)==3 || H(i,3)==4 ||  H(i,3)==5 || H(i,3)==6)
        H1stOrder(i,:)=H(i,:);
    end
end
H1stOrder(all(~H1stOrder,2),:)=[]; %Eliminate zero rows 
H1stOrder=sortrows(H1stOrder,-5);

for i=1:L(1)
    if abs(H(i,2))<SFTol || abs(H(i,4))<SFTol   %1st Filter condition:ignores rows with SF<0.3 
       H(i,:)=zeros(1,6);
    end
end

for i=1:L(1)
    if abs(H(i,6))<0.6    %2nd Filter condition:ignores rows with abs(m-prime)<0.3
        H(i,:)=zeros(1,6);
    end
end

H(all(~H,2),:)=[]; %Eliminate the zero rows
H=sortrows(H,-6); %Sort matrix in descending order of mprime values
HP=[H,abs(H(:,end))];
HP=sortrows(HP,-7);
%{ 
%Suppressed Tabulation
%Tabulate the information
 f = figure('Position', [100 100 752 350]);
 t = uitable('Parent', f, 'Position', [75 75 700 300]);
 set(t,'Data',HP);
 set(t,'ColumnName', {'SSnumber_Grain1','SchmidFactor_Grain1','SSnumber_Grain2','SchmidFactor_Grain2','Avg_SF','m_prime','abs(mp)'});
 
%  disp(H);
%Tabulate the basal and 1st Order prismatic slip systems:essentially a
%sanity check
j=figure('Position',[100 100 752 350]);
rt=uitable('Parent',j,'Position',[75 75 700 300]);
set(rt,'Data',H1stOrder);
set(rt,'ColumnName', {'SSnumber_Grain1','SchmidFactor_Grain1','SSnumber_Grain2','SchmidFactor_Grain2','Avg_SF','m_prime'});
%}



end