% Function to calculate average misorientation of two hcp crystals given
% Euler angles
function MisOr=HexMisOr(Euler1,Euler2)
g1=OrMat(Euler1);
g2=OrMat(Euler2);
g=g2/g1;
GM=HexSymm(g);
Ang=zeros(1,12); % Initialization of Angle Array

 % Rounding to 2 Decimal places is done in order to ensure that values stay within [-1,1] interval!
    Ang(1)=acosd((GM(1,1,1)+GM(2,2,1)+GM(3,3,1)-1)*0.5);
    Ang(2)=acosd((GM(1,1,2)+GM(2,2,2)+GM(3,3,2)-1)*0.5);
    Ang(3)=acosd((GM(1,1,3)+GM(2,2,3)+GM(3,3,3)-1)*0.5);
    Ang(4)=acosd((GM(1,1,4)+GM(2,2,4)+GM(3,3,4)-1)*0.5);
    Ang(5)=acosd((GM(1,1,5)+GM(2,2,5)+GM(3,3,5)-1)*0.5);
    Ang(6)=acosd((GM(1,1,6)+GM(2,2,6)+GM(3,3,6)-1)*0.5);
    Ang(7)=acosd((GM(1,1,7)+GM(2,2,7)+GM(3,3,7)-1)*0.5);
    Ang(8)=acosd((GM(1,1,8)+GM(2,2,8)+GM(3,3,8)-1)*0.5);
    Ang(9)=acosd((GM(1,1,9)+GM(2,2,9)+GM(3,3,9)-1)*0.5);
    Ang(10)=acosd((GM(1,1,10)+GM(2,2,10)+GM(3,3,10)-1)*0.5);
    Ang(11)=acosd((GM(1,1,11)+GM(2,2,11)+GM(3,3,11)-1)*0.5);
    Ang(12)=acosd((GM(1,1,12)+GM(2,2,12)+GM(3,3,12)-1)*0.5);

    MisOr=min(abs(Ang));
      





end



 
   
    
  
