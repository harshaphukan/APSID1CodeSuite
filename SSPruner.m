%Function to analyze and shortlist compatible slip system w.r.t a twin
% parent grain
function [F]=SSPruner(A,EulerP,STensP,ActiveVar,NIndex)
% Inputs
%A-->Array with Grain#,COM,Orientation,StressTensor for layer 
% LA-->Nearest neighbor grain information in same or different layer
% EulerP-->Orientation of parent grain
% ActiveVar-->Index corresponding to activated twin variant (between 34-39)
%TwinCOM--> Coordinates of the Twin COM
%NIndex--> Index of Neighbor Grain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine  Slip/Twin Systems that are within 70degrees w.r.t Vector connecting 
%COM of Twin and Neighbor in concern! Takes One Layer @ a time!
% Same Layer First!
NSDL=zeros(57,3); % Stacked array to store Slip/Twin  vectors for each neighbor insame layer 
SSList=zeros(1,57); % Array to store Slip/Twin System # that meets criteria
[sortmvL,slpsysL]=HexEval(A(NIndex,10:12),A(NIndex,4:9));
SNL=SlipSysFun(sortmvL,slpsysL);
for j=1:57
NSDL(j,:)=SNL(j).f3;
%if (acosd(dot(LATwin(i,2:4),NSDSL(j,:,i)))>=0) && (acosd(dot(LATwin(i,2:4),NSDSL(j,:,i)))<=70)
SSList(1,j)=SNL(j).f1;


%end

end
[H]=MPrimeFunHJP(EulerP,A(NIndex,10:12),STensP,A(NIndex,4:9),0.2,1e-6);



for i=1:length(H)
if  (H(i,end)>0.70) && (H(i,1)==ActiveVar)
    F(i,:)=H(i,:);
else
    F(i,:)=zeros(1,7);
    
end
end
F= F(any(F,2),:);

    

end


