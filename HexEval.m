%% Function to evaluate Schmid Factor of slip/twin systems, and respective
%% plane normals,directions
%% works only for hexagonal system [Based on code developed by Dr. T. R Bieler]
function [sortmv,slpsys]=HexEval(Euler,STens)
%Input1: EulerAngle (tuple of phi1,phi,phi2)
%StressTensor in Sample Coordinate system (Sxx,Syy,Szz,Sxy,Sxz,Syz) tuple
Ang=[1, 0, Euler(1),Euler(2),Euler(3)];

sigma=[STens(1),STens(4),STens(5);
       STens(4),STens(2),STens(6);
       STens(5),STens(6),STens(3)];% Initialize array to store stress tensors in sample coordinate system


str2 = sigma*sigma';
sigma_mag = (str2(1,1)+str2(2,2)+str2(3,3))^.5;
sigma_n= sigma/sigma_mag;   % normalized stress tensor to get generalized Schmid factor

c_a = 1.587;   nslphex = 57;    slpsys = cell(4,nslphex);    
%  Below, four points define plane, start point is start of b, 3rd is end
%  1st to 2nd or 2nd to 3rd or 3rd to 4th cross product identifines plane normal 
% Mode 1,  plane; direction, Define four points in the plane 

O = [ 0  0  0 0];   %          (I)         (H)
A = [ 2 -1 -1 0]/3; %            C ------- B 
B = [ 1  1 -2 0]/3; %          /  \       / \
C = [-1  2 -1 0]/3; %         /    a2   /    \
D = [-2  1  1 0]/3; %        /       \ /      \
E = [-1 -1  2 0]/3; %    (J)D ------ O(P)-a1-> A(G) 
F = [ 1 -2  1 0]/3; %        \       / \      /
P = [ 0  0  0 1];   %         \    a3   \    /
G = [ 2 -1 -1 3]/3; %          \  /       \ /
H = [ 1  1 -2 3]/3; %           E -------- F
I = [-1  2 -1 3]/3; %         (K)          (L)
J = [-2  1  1 3]/3; %
K = [-1 -1  2 3]/3; % where a1 = OA  a2 = OC  a3 = OE
L = [ 1 -2  1 3]/3; % DA || a1,  FC || a2,  BE || a3

% Slip system definitions set as of 1 July 2013 to be consistent with DAMASK
%   Hex           plane   direction   1st and 4th point is Burgers vector
%  basal <a>-glide:&Be Mg Re Ti; Re
slpsys{2,1} = [0 0 0 1;  2 -1 -1 0; D ; E ; F ; A ; B ; C ]; 
slpsys{2,2} = [0 0 0 1; -1  2 -1 0; F ; A ; B ; C ; D ; E ]; 
slpsys{2,3} = [0 0 0 1; -1 -1  2 0; B ; C ; D ; E ; F ; A ]; 

%  prism <a>-glide:Ti Zr RE; Be Re Mg
slpsys{2,4} = [ 0  1 -1 0;  2 -1 -1 0; E ; E ; E ; F ; L ; K ];
slpsys{2,5} = [-1  0  1 0; -1  2 -1 0; A ; A ; A ; B ; H ; G ];
slpsys{2,6} = [ 1 -1  0 0; -1 -1  2 0; C ; C ; C ; D ; J ; I ];

% prism <aa>:2nd order prismatic
slpsys{2,7} = [ 2 -1 -1 0;  0  1 -1 0; F; F; F; B ; H ; L ];
slpsys{2,8} = [-1  2 -1 0; -1  0  1 0; B; B; B; D ; J ; H ];
slpsys{2,9} = [-1 -1  2 0;  1 -1  0 0; D; D; D; F ; L ; J ];

%  pyramidal <a>-glide
slpsys{2,10} = [ 0 -1  1 1;  2 -1 -1 0; E ; E ; E ; F ; G ; J ];
slpsys{2,11} = [ 1  0 -1 1; -1  2 -1 0; A ; A ; A ; B ; I ; L ];
slpsys{2,12} = [-1  1  0 1; -1 -1  2 0; C ; C ; C ; D ; K ; H ];
slpsys{2,13} = [ 1 -1  0 1;  1  1 -2 0; F ; F ; F ; A ; H ; K ];
slpsys{2,14} = [ 0  1 -1 1; -2  1  1 0; B ; B ; B ; C ; J ; G ];
slpsys{2,15} = [-1  0  1 1;  1 -2  1 0; D ; D ; D ; E ; L ; I ];

%  pyramidal <c+a>-glide:; all?
slpsys{2,16} = [-1  1  0 1;  2 -1 -1 3; D ; D ; K ; P ; H ; C ]; 
slpsys{2,17} = [-1  1  0 1;  1 -2  1 3; C ; D ; K ; P ; H ; C ];
slpsys{2,18} = [ 1  0 -1 1; -1 -1  2 3; B ; B ; I ; P ; L ; A ];
slpsys{2,19} = [ 1  0 -1 1; -2  1  1 3; A ; B ; I ; P ; L ; A ];
slpsys{2,20} = [ 0 -1  1 1; -1  2 -1 3; F ; F ; G ; P ; J ; E ];
slpsys{2,21} = [ 0 -1  1 1;  1  1 -2 3; E ; F ; G ; P ; J ; E ];
slpsys{2,22} = [ 1 -1  0 1; -2  1  1 3; A ; A ; H ; P ; K ; F ];
slpsys{2,23} = [ 1 -1  0 1; -1  2 -1 3; F ; A ; H ; P ; K ; F ];
slpsys{2,24} = [-1  0  1 1;  1  1 -2 3; E ; E ; L ; P ; I ; D ];
slpsys{2,25} = [-1  0  1 1;  2 -1 -1 3; D ; E ; L ; P ; I ; D ]; 
slpsys{2,26} = [ 0  1 -1 1;  1 -2  1 3; C ; C ; J ; P ; G ; B ];
slpsys{2,27} = [ 0  1 -1 1; -1 -1  2 3; B ; C ; J ; P ; G ; B ];

%  pyramidal <c+a>-2nd order glide
slpsys{2,28} = [-2  1  1 1;  2 -1 -1 3; (O+D)/2 ; E ; L ; (P+G)/2 ; H ; C]; 
slpsys{2,29} = [ 1 -2  1 1; -1  2 -1 3; (O+F)/2 ; A ; H ; (P+I)/2 ; J ; E]; 
slpsys{2,30} = [ 1  1 -2 1; -1 -1  2 3; (O+B)/2 ; C ; J ; (P+K)/2 ; L ; A]; 
slpsys{2,31} = [ 2 -1 -1 1; -2  1  1 3; (O+A)/2 ; B ; I ; (P+J)/2 ; K ; F]; 
slpsys{2,32} = [-1  2 -1 1;  1 -2  1 3; (O+C)/2 ; D ; K ; (P+L)/2 ; G ; B]; 
slpsys{2,33} = [-1 -1  2 1;  1  1 -2 3; (O+E)/2 ; F ; G ; (P+H)/2 ; I ; D]; 

% *** Twin directions are opposite in Christian and Mahajan, and are not correcte to to be consistent with them
% FROM Kocks SXHEX   plane   direction   1st and 4th point is Burgers vector, order of C1 differs from Kock's file
%  {1012}<1011> T1 twins 0.17; -1.3  twins: all  Twin Vector must go in the
%  sense of shear, opposite C&M sense.
slpsys{2,34} = [-1  1  0 2;  1 -1  0 1; C ; D ; L ; G ; G ; G ]; 
slpsys{2,35} = [ 1  0 -1 2; -1  0  1 1; A ; B ; J ; K ; K ; K ]; 
slpsys{2,36} = [ 0 -1  1 2;  0  1 -1 1; E ; F ; H ; I ; I ; I ]; 
slpsys{2,37} = [ 1 -1  0 2; -1  1  0 1; F ; A ; I ; J ; J ; J ]; 
slpsys{2,38} = [-1  0  1 2;  1  0 -1 1; D ; E ; G ; H ; H ; H ]; 
slpsys{2,39} = [ 0  1 -1 2;  0 -1  1 1; B ; C ; K ; L ; L ; L ]; 

%  {2111}<2116> T2 twins: 0.63;  -0.4; Ti Zr Re RE]; Also does not follow C&M definition for shear direction
slpsys{2,40} = [-2  1  1 1;  2 -1 -1 6; (O+D)/2 ; E ; (L+K)/2 ; P ; (H+I)/2 ; C]; 
slpsys{2,41} = [ 1 -2  1 1; -1  2 -1 6; (O+F)/2 ; A ; (H+G)/2 ; P ; (J+K)/2 ; E]; 
slpsys{2,42} = [ 1  1 -2 1; -1 -1  2 6; (O+B)/2 ; C ; (J+I)/2 ; P ; (L+G)/2 ; A]; 
slpsys{2,43} = [ 2 -1 -1 1; -2  1  1 6; (O+A)/2 ; B ; (I+H)/2 ; P ; (K+L)/2 ; F]; 
slpsys{2,44} = [-1  2 -1 1;  1 -2  1 6; (O+C)/2 ; D ; (K+J)/2 ; P ; (G+H)/2 ; B]; 
slpsys{2,45} = [-1 -1  2 1;  1  1 -2 6; (O+E)/2 ; F ; (G+L)/2 ; P ; (I+J)/2 ; D]; 

%   {1011}<101-2> C1 twins: 0.10; 1.1; Mg; Zr Ti]; agrees with C&M
slpsys{2,46} = [-1  1  0 1; -1  1  0 -2; P ; H ; C ; (C+D)/2 ; D ; K ]; 
slpsys{2,47} = [ 1  0 -1 1;  1  0 -1 -2; P ; L ; A ; (A+B)/2 ; B ; I ]; 
slpsys{2,48} = [ 0 -1  1 1;  0 -1  1 -2; P ; J ; E ; (E+F)/2 ; F ; G ]; 
slpsys{2,49} = [ 1 -1  0 1;  1 -1  0 -2; P ; K ; F ; (F+A)/2 ; A ; H ]; 
slpsys{2,50} = [-1  0  1 1; -1  0  1 -2; P ; I ; D ; (D+E)/2 ; E ; L ]; 
slpsys{2,51} = [ 0  1 -1 1;  0  1 -1 -2; P ; G ; B ; (B+C)/2 ; C ; J ]; 

%  {2112}<211-3> C2 twins:; 0.22; 1.2 Ti Zr Re]; agrees with C&M
slpsys{2,52} = [ 2 -1 -1 2;  2 -1 -1 -3; (J+P)/2 ; K ; F ; (F+B)/2 ; B ; I]; 
slpsys{2,53} = [-1  2 -1 2; -1  2 -1 -3; (L+P)/2 ; G ; B ; (B+D)/2 ; D ; K]; 
slpsys{2,54} = [-1 -1  2 2; -1 -1  2 -3; (H+P)/2 ; I ; D ; (D+F)/2 ; F ; G]; 
slpsys{2,55} = [-2  1  1 2; -2  1  1 -3; (G+P)/2 ; H ; C ; (C+E)/2 ; E ; L]; 
slpsys{2,56} = [ 1 -2  1 2;  1 -2  1 -3; (I+P)/2 ; J ; E ; (E+A)/2 ; A ; H]; 
slpsys{2,57} = [ 1  1 -2 2;  1  1 -2 -3; (K+P)/2 ; L ; A ; (A+C)/2 ; C ; J]; 



nslphex=57;  %f2pyrc
c_a=1.587;

%
for isc=1:1:nslphex  % cell 3 = slip system unit vectors, cell 4 = Schmid matrix
    n=[slpsys{2,isc}(1,1), (slpsys{2,isc}(1,1)+2*slpsys{2,isc}(1,2))/sqrt(3), slpsys{2,isc}(1,4)/c_a]; %plane normal in cartesian
    m=[3*slpsys{2,isc}(2,1)/2, (slpsys{2,isc}(2,1)+2*slpsys{2,isc}(2,2))*sqrt(3)/2, slpsys{2,isc}(2,4)*c_a];% slip direction in cartesian
    p1=[3*slpsys{2,isc}(3,1)/2, (slpsys{2,isc}(3,1)+2*slpsys{2,isc}(3,2))*sqrt(3)/2, slpsys{2,isc}(3,4)*c_a]; % point 1
    p2=[3*slpsys{2,isc}(4,1)/2, (slpsys{2,isc}(4,1)+2*slpsys{2,isc}(4,2))*sqrt(3)/2, slpsys{2,isc}(4,4)*c_a]; % point 2
    p3=[3*slpsys{2,isc}(5,1)/2, (slpsys{2,isc}(5,1)+2*slpsys{2,isc}(5,2))*sqrt(3)/2, slpsys{2,isc}(5,4)*c_a]; % point 3
    p4=[3*slpsys{2,isc}(6,1)/2, (slpsys{2,isc}(6,1)+2*slpsys{2,isc}(6,2))*sqrt(3)/2, slpsys{2,isc}(6,4)*c_a]; % point 4
    p5=[3*slpsys{2,isc}(7,1)/2, (slpsys{2,isc}(7,1)+2*slpsys{2,isc}(7,2))*sqrt(3)/2, slpsys{2,isc}(7,4)*c_a]; % point 5
    p6=[3*slpsys{2,isc}(8,1)/2, (slpsys{2,isc}(8,1)+2*slpsys{2,isc}(8,2))*sqrt(3)/2, slpsys{2,isc}(8,4)*c_a]; % point 6
    mag_m=(m(1,1)^2+m(1,2)^2+m(1,3)^2)^0.5; 
    mag_n=(n(1,1)^2+n(1,2)^2+n(1,3)^2)^0.5; 
    unit_m = m/mag_m;
    unit_n = n/mag_n;
    slpsys{1,isc} = isc;
    slpsys{3,isc} = [unit_n;unit_m;p1;p2;p3;p4;p5;p6];  % normal in first row, direction in next row, points in next 6 rows
    slpsys{4,isc} = 0.5*(unit_m'*unit_n + unit_n'*unit_m); % Schmid matrix
end

sij = [0.9581 -0.4623 -0.1893 0.698 2.1413 2.8408]/100; % for Ti from Simmons and Wang  in units of 1/GPa
Elow = 83.2;   Ehigh = 145.5;   % Elastic Ccontants of Nb from Simmons and Wang for Ti


lAng = size(Ang);

for iAng=1:1:lAng(1,1) %1;%4; %2; %  
    phidA = Ang(iAng,3:5);
    phidA(1) = phidA(1); %+180; % **** Rotating euler angles (e.g. to correct for 180 rotation)
    if phidA(1)>360   
        phidA(1) = phidA(1)-360; 
    end
    if phidA(1)<0            
        phidA(1) = phidA(1)+360; 
    end
    phisA = phidA*pi/180;  %Compute Bunge orientation matrix g
    gphi1=[cos(phisA(1,1)),sin(phisA(1,1)),0;-sin(phisA(1,1)),cos(phisA(1,1)),0;0,0,1]; 
    gPhi=[1,0,0;0,cos(phisA(1,2)),sin(phisA(1,2));0,-sin(phisA(1,2)),cos(phisA(1,2))];
    gphi2=[cos(phisA(1,3)),sin(phisA(1,3)),0;-sin(phisA(1,3)),cos(phisA(1,3)),0;0,0,1];
    gA=gphi2*gPhi*gphi1;    
    g(:,:,iAng)=gA;
    gsgTA = gA*sigma_n(:,:,iAng)*gA'; %rotated stress tensor
    
            
    for isc=1:1:nslphex %  Compute Schmid Factors, followed by grain slip system; slpsys{4 = Schmid matrix}
        computeSFA(isc,1)=isc;
        computeSFA(isc,2)=0.;
        for i=1:1:3
            for j=1:1:3
                computeSFA(isc,2)=computeSFA(isc,2)+gsgTA(i,j)*slpsys{4,isc}(i,j);
                
            end
        end               % Abs(Schmid factor) is in column after signed Schmid factor
        if isc>24 && computeSFA(isc,2)<0 
            computeSFA(isc,2) = 0.001*computeSFA(isc,2) ;    % this is to prevent anti-twin shears from being seriously considered later
        end
        computeSFA(isc,3:11)=[abs(computeSFA(isc,2)),slpsys{2,isc}(1,:),slpsys{2,isc}(2,:)];
        rot_nA = slpsys{3,isc}(1,:)*gA;
        rot_bA = slpsys{3,isc}(2,:)*gA;
        rot_p1 = slpsys{3,isc}(3,:)*gA;
        rot_p2 = slpsys{3,isc}(4,:)*gA;
        rot_p3 = slpsys{3,isc}(5,:)*gA;
        rot_p4 = slpsys{3,isc}(6,:)*gA;
        rot_p5 = slpsys{3,isc}(7,:)*gA;
        rot_p6 = slpsys{3,isc}(8,:)*gA;
% Schmid factors (1-3), rotated plane normal (4-6), Computed rotated Burgers vector (7-9),  
% , plane trace on Z surface (10-12) Computed rotated position vectors to points p1-p4 (13-24)
        Schm_labvecA(isc,:) = [computeSFA(isc,1:3), rot_nA, rot_bA, cross(rot_nA',[0,0,1]), rot_p1, rot_p2, rot_p3, rot_p4, rot_p5, rot_p6];
    end                             %    1-3         4-6      7-9               10-12        13-15   16-18   19-21   22-24   25-27   28-30
    %useful plotting unit cell vectors that sort to bottom row
    Schm_labvecA(nslphex+1,:) = [0 1 -1 [1 0 0]*gA [0 1 0]*gA [0 0 1]*gA 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
    sortmv(:,:,iAng) = sortrows(Schm_labvecA,-3);
end




end