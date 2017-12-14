% Function to calculate unit vectors linking a given point and a list of
% points(or one single pt) in 2D
function U=UnitVec(VArray,V)
if numel(VArray)==2
    U(1)=V(1)-VArray(1);
    U(2)=V(2)-VArray(2);
    U=U/norm(U);
else
    U=zeros(length(VArray),2);
for i=1:length(U)
      U(i,1)=V(1)-VArray(i,1);
      U(i,2)=V(2)-VArray(i,2);
      U(i,:)=U(i,:)/norm(U(i,:));
end
end

end