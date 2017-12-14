function g=OrMat(Euler)
% Function to calculate Orientation Matrix Given Euler Angles (degrees) as a vector
sphi1=sind(Euler(1));sphi=sind(Euler(2)); sphi2=sind(Euler(3));
cphi1=cosd(Euler(1));cphi=cosd(Euler(2)); cphi2=cosd(Euler(3));

gphi1=[cphi1,sphi1,0;-sphi1,cphi1,0;0,0,1];
gphi=[1,0,0;0,cphi,sphi;0,-sphi,cphi];
gphi2=[cphi2,sphi2,0;-sphi2,cphi2,0;0,0,1];

g=gphi2*gphi*gphi1;

       
end