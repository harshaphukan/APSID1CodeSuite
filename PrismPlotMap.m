%% Plots the hexagonal unit cells on the Center of Mass Positions on the GrainMap
function PrismPlotMap(Array)
%ptpl = 0;  % 0 = don't plot plane traces
[R,~]=size(Array);
N_grain=R;
%%
%% Part 2
for i=1:R
    Ang(i,:)=[i, 0, Array(i,10),Array(i,11),Array(i,12)];% Note Ang(1,:) corresponds to grain 1 and so on;
end

sigma=zeros(3,3,N_grain); % Initialize array to store stress tensors in sample coordinate system
for j=1:R
    
sigma(:,:,j)= [Array(j,4),Array(j,7),Array(j,8);Array(j,7),Array(j,5),Array(j,9);Array(j,8),Array(j,9),Array(j,6)]; 
end


nsten =N_grain;

for i = 1:1:nsten
    str2 = sigma(:,:,i)*sigma(:,:,i)';
    sigma_mag = (str2(1,1)+str2(2,2)+str2(3,3))^.5;
    sigma_n(:,:,i) = sigma(:,:,i)/sigma_mag;   % normalized stress tensor to get generalized Schmid factor
end
stereo = 1;  % plots stereographic projection, otherwise direct projection
ipfd = 2;  %inverse pole figure direction (x,y,z, = 1,2,3) will plot points from a group
% CTEd = ipfd; % plot magnitude of CTE in this direction as gray scale within symbol
iEd = ipfd;  % plot magnitude of E in this direction as red (high) - blue (low)
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
ibas = 1;
fbas = 3;
%  prism <a>-glide:Ti Zr RE; Be Re Mg
slpsys{2,4} = [ 0  1 -1 0;  2 -1 -1 0; E ; E ; E ; F ; L ; K ];
slpsys{2,5} = [-1  0  1 0; -1  2 -1 0; A ; A ; A ; B ; H ; G ];
slpsys{2,6} = [ 1 -1  0 0; -1 -1  2 0; C ; C ; C ; D ; J ; I ];
iprs = 4;
fprs = 6;
% prism <aa>:2nd order prismatic
slpsys{2,7} = [ 2 -1 -1 0;  0  1 -1 0; F; F; F; B ; H ; L ];
slpsys{2,8} = [-1  2 -1 0; -1  0  1 0; B; B; B; D ; J ; H ];
slpsys{2,9} = [-1 -1  2 0;  1 -1  0 0; D; D; D; F ; L ; J ];
i2prs = 7;
f2prs = 9;
%  pyramidal <a>-glide
slpsys{2,10} = [ 0 -1  1 1;  2 -1 -1 0; E ; E ; E ; F ; G ; J ];
slpsys{2,11} = [ 1  0 -1 1; -1  2 -1 0; A ; A ; A ; B ; I ; L ];
slpsys{2,12} = [-1  1  0 1; -1 -1  2 0; C ; C ; C ; D ; K ; H ];
slpsys{2,13} = [ 1 -1  0 1;  1  1 -2 0; F ; F ; F ; A ; H ; K ];
slpsys{2,14} = [ 0  1 -1 1; -2  1  1 0; B ; B ; B ; C ; J ; G ];
slpsys{2,15} = [-1  0  1 1;  1 -2  1 0; D ; D ; D ; E ; L ; I ];
ipyra = 10;
fpyra = 15;
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
ipyrc = 16;
fpyrc = 27;
%  pyramidal <c+a>-2nd order glide
slpsys{2,28} = [-2  1  1 1;  2 -1 -1 3; (O+D)/2 ; E ; L ; (P+G)/2 ; H ; C]; 
slpsys{2,29} = [ 1 -2  1 1; -1  2 -1 3; (O+F)/2 ; A ; H ; (P+I)/2 ; J ; E]; 
slpsys{2,30} = [ 1  1 -2 1; -1 -1  2 3; (O+B)/2 ; C ; J ; (P+K)/2 ; L ; A]; 
slpsys{2,31} = [ 2 -1 -1 1; -2  1  1 3; (O+A)/2 ; B ; I ; (P+J)/2 ; K ; F]; 
slpsys{2,32} = [-1  2 -1 1;  1 -2  1 3; (O+C)/2 ; D ; K ; (P+L)/2 ; G ; B]; 
slpsys{2,33} = [-1 -1  2 1;  1  1 -2 3; (O+E)/2 ; F ; G ; (P+H)/2 ; I ; D]; 
i2pyrc = 28;
f2pyrc = 33;
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
iT1 = 34;
fT1 = 39;
%  {2111}<2116> T2 twins: 0.63;  -0.4; Ti Zr Re RE]; Also does not follow C&M definition for shear direction
slpsys{2,40} = [-2  1  1 1;  2 -1 -1 6; (O+D)/2 ; E ; (L+K)/2 ; P ; (H+I)/2 ; C]; 
slpsys{2,41} = [ 1 -2  1 1; -1  2 -1 6; (O+F)/2 ; A ; (H+G)/2 ; P ; (J+K)/2 ; E]; 
slpsys{2,42} = [ 1  1 -2 1; -1 -1  2 6; (O+B)/2 ; C ; (J+I)/2 ; P ; (L+G)/2 ; A]; 
slpsys{2,43} = [ 2 -1 -1 1; -2  1  1 6; (O+A)/2 ; B ; (I+H)/2 ; P ; (K+L)/2 ; F]; 
slpsys{2,44} = [-1  2 -1 1;  1 -2  1 6; (O+C)/2 ; D ; (K+J)/2 ; P ; (G+H)/2 ; B]; 
slpsys{2,45} = [-1 -1  2 1;  1  1 -2 6; (O+E)/2 ; F ; (G+L)/2 ; P ; (I+J)/2 ; D]; 
iT2 = 40;
fT2 = 45;
%   {1011}<101-2> C1 twins: 0.10; 1.1; Mg; Zr Ti]; agrees with C&M
slpsys{2,46} = [-1  1  0 1; -1  1  0 -2; P ; H ; C ; (C+D)/2 ; D ; K ]; 
slpsys{2,47} = [ 1  0 -1 1;  1  0 -1 -2; P ; L ; A ; (A+B)/2 ; B ; I ]; 
slpsys{2,48} = [ 0 -1  1 1;  0 -1  1 -2; P ; J ; E ; (E+F)/2 ; F ; G ]; 
slpsys{2,49} = [ 1 -1  0 1;  1 -1  0 -2; P ; K ; F ; (F+A)/2 ; A ; H ]; 
slpsys{2,50} = [-1  0  1 1; -1  0  1 -2; P ; I ; D ; (D+E)/2 ; E ; L ]; 
slpsys{2,51} = [ 0  1 -1 1;  0  1 -1 -2; P ; G ; B ; (B+C)/2 ; C ; J ]; 
iC1 = 46;
fC1 = 51;
%  {2112}<211-3> C2 twins:; 0.22; 1.2 Ti Zr Re]; agrees with C&M
slpsys{2,52} = [ 2 -1 -1 2;  2 -1 -1 -3; (J+P)/2 ; K ; F ; (F+B)/2 ; B ; I]; 
slpsys{2,53} = [-1  2 -1 2; -1  2 -1 -3; (L+P)/2 ; G ; B ; (B+D)/2 ; D ; K]; 
slpsys{2,54} = [-1 -1  2 2; -1 -1  2 -3; (H+P)/2 ; I ; D ; (D+F)/2 ; F ; G]; 
slpsys{2,55} = [-2  1  1 2; -2  1  1 -3; (G+P)/2 ; H ; C ; (C+E)/2 ; E ; L]; 
slpsys{2,56} = [ 1 -2  1 2;  1 -2  1 -3; (I+P)/2 ; J ; E ; (E+A)/2 ; A ; H]; 
slpsys{2,57} = [ 1  1 -2 2;  1  1 -2 -3; (K+P)/2 ; L ; A ; (A+C)/2 ; C ; J]; 
iC2 = 52;
fC2 = 57;


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


%%
%  Strategy:  First extract useful vectors from slip system information to draw the hexagonal prisms
%  positions in sortmv  p1:13-15  p2:16-18  p3:19-21  p4:22-24  p5:25-27  p6:28-30 
%  positions in pln    p1:4-6   p2:7-9   p3:10-12  p4:13-15  p5:16-18  p6:18-21
plnM=zeros(8,21,N_grain);sortplnM=zeros(8,21,N_grain);
rotcM=zeros(N_grain,3);a1=zeros(N_grain,3);a2=zeros(N_grain,3);a3=zeros(N_grain,3); center=zeros(N_grain,3);
minxM=zeros(N_grain,1);minyM=zeros(N_grain,1); minzM=zeros(N_grain,1); maxxM=zeros(N_grain,1);maxyM=zeros(N_grain,1); maxzM=zeros(N_grain,1);
midxM=zeros(N_grain,1); midyM=zeros(N_grain,1);
prsxpltM=zeros(8,7,N_grain); prsypltM=zeros(8,7,N_grain); prszpltM=zeros(8,7,N_grain);
sp1M=zeros(N_grain,3);sp2M=zeros(N_grain,3);sp3M=zeros(N_grain,3);sp4M=zeros(N_grain,3);sp5M=zeros(N_grain,3);sp6M=zeros(N_grain,3); 
for GN=1:1:N_grain
for isc = 1:1:nslphex   
    if sortmv(isc,1,GN) == 1                % locate the two basal planes
        plnM(1,4:21,GN) = sortmv(isc,13:30,GN);  % bottom basal plane
        plnM(2,4:21,GN) = sortmv(isc,13:30,GN);  % top basal plane
        rotcM(GN,:) = sortmv(isc,4:6,GN)*c_a;       % basal plane normal * c/a
        for j = 4:3:19
            plnM(2,j:j+2,GN) = plnM(1,j:j+2,GN) + rotcM(GN,:);    % move top plane up by a unit of c
        end
        a1(GN,:)= sortmv(isc,7:9,GN);   %  locate a1 using SS1
    elseif sortmv(isc,1,GN) == 2
        a2(GN,:) = sortmv(isc,7:9,GN);   %  locate a2 using SS2
    elseif sortmv(isc,1,GN) == 3
        a3(GN,:) = sortmv(isc,7:9,GN);   %  locate a3 using SS3
    end
end
for isc = 1:1:nslphex
    n = [0 0 0 sortmv(isc,4:6,iAng)];
    if sortmv(isc,1,GN) == 4       %  locate two prism planes on opposite sides using SS4
    	plnM(3,4:21,GN) = sortmv(isc,13:30,GN);
        for j = 13:3:28
            plnM(4,j-9:j-7,GN) = sortmv(isc,j:j+2,GN) + a2(GN,:) - a3(GN,:);
        end
    elseif sortmv(isc,1,GN) == 5    %  locate two prism planes on opposite sides using SS5
    	plnM(5,4:21,GN) = sortmv(isc,13:30,GN);
     	for j = 13:3:28
            plnM(6,j-9:j-7,GN) = sortmv(isc,j:j+2,GN) + a3(GN,:) - a1(GN,:);
        end
    elseif sortmv(isc,1,GN) == 6   %  locate two prism planes on opposite sides using SS6
    	plnM(7,4:21,GN) = sortmv(isc,13:30,GN);
     	for j = 13:3:28
            plnM(8,j-9:j-7,GN) = sortmv(isc,j:j+2,GN) + a1(GN,:) - a2(GN,:);
        end
    end
end

for j = 1:1:2    % Find z elevation of basal planes
    for k = 1:1:3
        plnM(j,k,GN) = (plnM(j,3+k,GN)+plnM(j,6+k,GN)+plnM(j,9+k,GN)+plnM(j,12+k,GN)+plnM(j,15+k,GN)+plnM(j,18+k,GN))/6;
    end
end
center(GN,:) = (plnM(1,1:3,GN)+plnM(2,1:3,GN))/2; 
for j = 3:1:8  % Find z elevation of prism planes
    plnM(j,3,GN) = (plnM(j,12,GN)+plnM(j,15,GN)+plnM(j,18,GN)+plnM(j,21,GN))/4;
end
sortplnM(:,:,GN) = sortrows(plnM(:,:,GN),-3);
%minx = 0; miny = 0; minz = 0; maxx = 0; maxy = 0; maxz = 0;
for j = 1:1:8 % assemble vectors for plotting faces of hex prism
    prsxpltM(j,1:7,GN) = [sortplnM(j,4,GN) sortplnM(j,7,GN) sortplnM(j,10,GN) sortplnM(j,13,GN) sortplnM(j,16,GN) sortplnM(j,19,GN) sortplnM(j,4,GN)];
    minxM(GN) = min(minxM(GN),min(prsxpltM(j,:,GN))); maxxM(GN)= max(maxxM(GN),max(prsxpltM(j,:,GN)));
    prsypltM(j,1:7,GN) = [sortplnM(j,5,GN) sortplnM(j,8,GN) sortplnM(j,11,GN) sortplnM(j,14,GN) sortplnM(j,17,GN) sortplnM(j,20,GN) sortplnM(j,5,GN)];
    minyM(GN) = min(minyM(GN),min(prsypltM(j,:,GN))); maxyM(GN) = max(maxyM(GN),max(prsypltM(j,:,GN)));
    prszpltM(j,1:7,GN) = [sortplnM(j,6,GN) sortplnM(j,9,GN) sortplnM(j,12,GN) sortplnM(j,15,GN) sortplnM(j,18,GN) sortplnM(j,21,GN) sortplnM(j,6,GN)];
    minzM(GN)= min(minzM(GN),min(prszpltM(j,:,GN))); maxzM(GN) = max(maxzM(GN),max(prszpltM(j,:,GN)));
end
    sp1M(GN,:) = sortmv(1,13:15,GN); % identify plotted points on the slip plane
    sp2M(GN,:) = sortmv(1,16:18,GN);
    sp3M(GN,:) = sortmv(1,19:21,GN);
    sp4M(GN,:) = sortmv(1,22:24,GN);
    sp5M(GN,:) = sortmv(1,25:27,GN);
    sp6M(GN,:) = sortmv(1,28:30,GN);
    minxM(GN) = min(minxM(GN), n(4));
    minyM(GN) = min(minyM(GN), n(5));
    maxxM(GN) = max(maxxM(GN), n(4));   % find appropriate range of x and y for plot
    maxyM(GN) = max(maxyM(GN), n(5));
    midxM(GN) = (minxM(GN)+maxxM(GN))/2;
    midyM(GN)= (minyM(GN)+maxyM(GN))/2;
    del = 1.4;
end
%%
%% This part rotates the COM Map by 25 degrees (CCW)
XM=zeros(N_grain,1);
YM=zeros(N_grain,1);
XY25=zeros(N_grain,2);
theta=deg2rad(25); % Rotate COM Positions by 30 degrees CCW in the XY Plane
%Define Orientation Matrix in 2D
RotMat=[cos(theta),-sin(theta);sin(theta),cos(theta)];

for i=1:N_grain
  AxRot=RotMat*[Array(i,2),Array(i,3)]';
  XY25(i,1)=AxRot(1);
  XY25(i,2)=AxRot(2);
end

%% The portion is to plot the spatial COMs of the grains

figure
axis square
axis('off')


title('Grain Map','Fontsize',25)
for count=1:N_grain 
    XM(count)=XY25(count,1)+800;%Rowstrn(count,7)+800;
    XM(count)=XM(count)/1400;
    YM(count)=XY25(count,2)+700;%Rowstrn(count,8)+700;
    YM(count)=YM(count)/1400;
    

axes('Position', [XM(count) YM(count) .05 .05], 'Layer','top');  %Normalized Position
    
    %% These plots will match TSL with X down

      axis square
%      axis([-1,1,-1,1],'off')
   set(gca ,'ycolor' ,'w'); set(gca ,'xcolor' ,'w');  % make axes white for ease in later arranging.  
   
    for j = 1:1:4 % plot the 4 top most surface prisms of the hex cell that have the highest z elevation
        
        hold on
        plot(prsypltM(j,:,count),-prsxpltM(j,:,count), 'Linewidth',2,'Color',[.0 .0 .0]);
        
    end
  
  text(XM(count)+0.5, YM(count),num2str(Array(count,1)),'Fontsize',10); 
 %text(XM(count)+0.5, YM(count),num2str(count),'Fontsize',10); 
 axis([midyM(count)-del midyM(count)+del -midxM(count)-del -midxM(count)+del],'off','tight') %Make axes invisible for the overlay
end
end

