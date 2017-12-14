% Function to calculate von Mises Stress for a Grain List in Layer
% Input: Array output by LoadData function for Beamline 1 Data set
function [vMArray]=vonMises(A)

L=length(A);
vMArray=zeros(L,1);
for i=1:L
    SigHydro=sum(A(i,1:6))/3.0;
    S11=A(i,1)-SigHydro; S22=A(i,2)-SigHydro;S33=A(i,3)-SigHydro;
    Sdev=[S11,A(i,4),A(i,5);A(i,4),S22,A(i,6);A(i,5),A(i,6),S33];
    vMArray(i)=sqrt(1.5*trace(Sdev*Sdev));
    
end
    
%stress(1,1),stress(2,2),stress(3,3),stress(1,2),stress(1,3),stress(2,3)



end