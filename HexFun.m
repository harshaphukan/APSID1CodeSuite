function HexFun(EulerAng,StressTens)
Ang(1,:) = [1 0 EulerAng(1) EulerAng(2) EulerAng(3)]; %38.7349551454330,63.2055056896633,198.434948822922   126, 85, 229 106 102 27  

sEul = size(Ang);
%
%  Stress tensor is defined; uncomment the one you want, or make a new one.  
%ststens = [1 0 0 ; 0 0 0 ; 0 0 0]; % tension in X  %  
%ststens = [0 0 0 ; 0 1 0 ; 0 0 0]; % tension in Y  %  
ststens = [StressTens(1), StressTens(4),StressTens(5); 
           StressTens(4),StressTens(2), StressTens(6); 
           StressTens(5),StressTens(6),StressTens(3)]; % tension in Z  %  
%ststens = [0 0 0 ; 0 0 -1 ; 0 -1 0]; % Shear in XZ plane  %  
str2 = ststens*ststens';
ststens_mag = (str2(1,1)+str2(2,2)+str2(3,3))^.5;
ststens_n = ststens/ststens_mag;   % normalized stress tensor to get generalized Schmid factor
%
%  Inverse pole figures can be drawn by setting values of ipfd; ipfd is
%  inverse pole figure direction, 
%  iCTE is CTE direction, iEd is E direction, not checked recently, so bugs may exist.
stereo = 1;  % plots stereographic projection, otherwise direct projection
ipfd = 3;  %inverse pole figure direction (x,y,z, = 1,2,3) will plot points from a group
ihex = 1;  % inverse pole figure for hexagonal crystal
iEd = ipfd;  % plot magnitude of E in this direction as red (high) - blue (low)
% CTEd = ipfd; % plot magnitude of CTE in this direction as gray scale within symbol
%
%  This section allows the point of view of the unit cell to be changed.
Rotateview = [1 0 0 ; 0 1 0 ; 0 0 1];   
iRotatev = 0;  % Provide NUMBER of rotations
% rotation of observer from normal TSL point of view to another point of view (crystal stays put)
% This right hand rotation matrix has a 4th row which has in columns 2,3, the angle and axis of the rotation.  
%Rotation(4,2:3,1) = [5,3];   % rotation to move axis to a tilted direction in the X-Y plane
Rotation(4,2:3,1) = [180,2];   % rotate viewpoint about vertical (x) axis + to view from above/right, - to view from above/left
%Rotation(4,2:3,2) = [90,2];   % +90 to view from below, -90 to view from above
% Build the Euler angle rotation matrix starting with identity
if iRotatev ~= 0
    for i = 1:1:iRotatev
    rang = Rotation(4,2,i);
    if Rotation(4,3,i) == 3
        Rotation(1:3,1:3,i)=[cosd(rang),sind(rang),0;-sind(rang),cosd(rang),0;0,0,1];
    elseif Rotation(4,3,i) == 2
        Rotation(1:3,1:3,i)=[cosd(rang),0,sind(rang);0,1,0;-sind(rang),0,cosd(rang)];
    elseif Rotation(4,3,i) == 1
        Rotation(1:3,1:3,i)=[1,0,0;0,cosd(rang),sind(rang);0,-sind(rang),cosd(rang)];
    end
    Rotateview = Rotateview*Rotation(1:3,1:3,i);
    end
end
%
%  Below, six points define plane, start point is start of Burgers vector b, 4th is end of b
%  1st to 2nd or 2nd to 3rd or 3rd to 4th cross product identifines plane normal 
nslphex=57; %69 
c_a=1.587;   % this is the last place for user input in this cell ...

slpsys = cell(4,nslphex);    

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

% Mode 1,  plane; direction, Define six points in the plane 
% basal <a>
slpsys{2,1} = [0 0 0 1;  2 -1 -1 0; D ; E ; F ; A ; B ; C ];
slpsys{2,2} = [0 0 0 1; -1  2 -1 0; F ; A ; B ; C ; D ; E ];
slpsys{2,3} = [0 0 0 1; -1 -1  2 0; B ; C ; D ; E ; F ; A ];

% prism <a>
slpsys{2,4} = [ 0 -1  1 0;  2 -1 -1 0; E ; E ; E ; F ; L ; K ];
slpsys{2,5} = [ 1  0 -1 0; -1  2 -1 0; A ; A ; A ; B ; H ; G ];
slpsys{2,6} = [ 1 -1  0 0; -1 -1  2 0; C ; C ; C ; D ; J ; I ];

% prism <aa>
slpsys{2,7} = [ 2 -1 -1 0; 0 -1  1 0; F; F; F; B ; H ; L ];
slpsys{2,8} = [-1  2 -1 0; 1  0 -1 0; B; B; B; D ; J ; H ];
slpsys{2,9} = [-1 -1  2 0; 1 -1  0 0; D; D; D; F ; L ; J ];

%pyramidal <a>
slpsys{2,10} =  [ 1  0 -1 1; -1  2 -1 0; A ; A ; A ; B ; I ; L ];
slpsys{2,11} =  [ 0  1 -1 1; -2  1  1 0; B ; B ; B ; C ; J ; G ];
slpsys{2,12} =  [-1  1  0 1; -1 -1  2 0; C ; C ; C ; D ; K ; H ];
slpsys{2,13} = [-1  0  1 1;  1 -2  1 0; D ; D ; D ; E ; L ; I ];
slpsys{2,14} = [ 0 -1  1 1;  2 -1 -1 0; E ; E ; E ; F ; G ; J ];
slpsys{2,15} = [ 1 -1  0 1;  1  1 -2 0; F ; F ; F ; A ; H ; K ];

%  pyramidal <c+a>-glide
slpsys{2,16} = [ 1  0 -1 1; -2  1  1 3; A ; B ; I ; P ; L ; A ];
slpsys{2,17} = [ 1  0 -1 1; -1 -1  2 3; B ; B ; I ; P ; L ; A ];
slpsys{2,18} = [ 0  1 -1 1; -1 -1  2 3; B ; C ; J ; P ; G ; B ];
slpsys{2,19} = [ 0  1 -1 1;  1 -2  1 3; C ; C ; J ; P ; G ; B ];
slpsys{2,20} = [-1  1  0 1;  1 -2  1 3; C ; D ; K ; P ; H ; C ];
slpsys{2,21} = [-1  1  0 1;  2 -1 -1 3; D ; D ; K ; P ; H ; C ]; 
slpsys{2,22} = [-1  0  1 1;  2 -1 -1 3; D ; E ; L ; P ; I ; D ]; 
slpsys{2,23} = [-1  0  1 1;  1  1 -2 3; E ; E ; L ; P ; I ; D ];
slpsys{2,24} = [ 0 -1  1 1;  1  1 -2 3; E ; F ; G ; P ; J ; E ];
slpsys{2,25} = [ 0 -1  1 1; -1  2 -1 3; F ; F ; G ; P ; J ; E ];
slpsys{2,26} = [ 1 -1  0 1; -1  2 -1 3; F ; A ; H ; P ; K ; F ];
slpsys{2,27} = [ 1 -1  0 1; -2  1  1 3; A ; A ; H ; P ; K ; F ];

%  pyramidal <c+a>-2nd order glide
slpsys{2,28} = [ 2 -1 -1 2; -2  1  1 3; (O+A)/2 ; B ; I ; (P+J)/2 ; K ; F]; 
slpsys{2,29} = [ 1  1 -2 2; -1 -1  2 3; (O+B)/2 ; C ; J ; (P+K)/2 ; L ; A]; 
slpsys{2,30} = [-1  2 -1 2;  1 -2  1 3; (O+C)/2 ; D ; K ; (P+L)/2 ; G ; B]; 
slpsys{2,31} = [-2  1  1 2;  2 -1 -1 3; (O+D)/2 ; E ; L ; (P+G)/2 ; H ; C]; 
slpsys{2,32} = [-1 -1  2 2;  1  1 -2 3; (O+E)/2 ; F ; G ; (P+H)/2 ; I ; D]; 
slpsys{2,33} = [ 1 -2  1 2; -1  2 -1 3; (O+F)/2 ; A ; H ; (P+I)/2 ; J ; E]; 

% FROM Kocks SXHEX   plane   direction   1st and 4th point is Burgers vector, order of C1 differs from Kock's file
%  {1012}<1011> T1 twins 0.17; -1.3  twins: all
slpsys{2,34} = [ 1  0 -1 2; -1  0  1 1; A ; B ; J ; K ; K ; K ]; 
slpsys{2,35} = [ 0  1 -1 2;  0 -1  1 1; B ; C ; K ; L ; L ; L ]; 
slpsys{2,36} = [-1  1  0 2;  1 -1  0 1; C ; D ; L ; G ; G ; G ]; 
slpsys{2,37} = [-1  0  1 2;  1  0 -1 1; D ; E ; G ; H ; H ; H ]; 
slpsys{2,38} = [ 0 -1  1 2;  0  1 -1 1; E ; F ; H ; I ; I ; I ]; 
slpsys{2,39} = [ 1 -1  0 2; -1  1  0 1; F ; A ; I ; J ; J ; J ]; 

%  <c+a>-3rd order glide; 
slpsys{2,40} = [ 2 -1 -1 1; -1  2 -1 3; F ; (K+L)/2 ; P ; P ; (I+H)/2 ; B ];  
slpsys{2,41} = [ 2 -1 -1 1; -1 -1  2 3; B ; (I+H)/2 ; P ; P ; (K+L)/2 ; F ];
slpsys{2,42} = [ 1  1 -2 1; -2  1  1 3; A ; (L+G)/2 ; P ; P ; (J+I)/2 ; C ]; 
slpsys{2,43} = [ 1  1 -2 1;  1 -2  1 3; C ; (J+I)/2 ; P ; P ; (L+G)/2 ; A ]; 
slpsys{2,44} = [-1  2 -1 1; -1 -1  2 3; B ; (G+H)/2 ; P ; P ; (K+J)/2 ; D ]; 
slpsys{2,45} = [-1  2 -1 1;  2 -1 -1 3; D ; (K+J)/2 ; P ; P ; (G+H)/2 ; B ]; 
slpsys{2,46} = [-2  1  1 1;  1 -2  1 3; C ; (H+I)/2 ; P ; P ; (L+K)/2 ; E ]; 
slpsys{2,47} = [-2  1  1 1;  1  1 -2 3; E ; (L+K)/2 ; P ; P ; (H+I)/2 ; C ]; 
slpsys{2,48} = [-1 -1  2 1;  2 -1 -1 3; D ; (I+J)/2 ; P ; P ; (G+L)/2 ; F ]; 
slpsys{2,49} = [-1 -1  2 1; -1  2 -1 3; F ; (G+L)/2 ; P ; P ; (I+J)/2 ; D ]; 
slpsys{2,50} = [ 1 -2  1 1;  1  1 -2 3; E ; (J+K)/2 ; P ; P ; (H+G)/2 ; A ]; 
slpsys{2,51} = [ 1 -2  1 1; -2  1  1 3; A ; (H+G)/2 ; P ; P ; (J+K)/2 ; E ]; 

%  {2111}<2116> T2 twins: 0.63;  -0.4; Ti Zr Re RE]; 
slpsys{2,52} = [ 2 -1 -1 1; -2  1  1 6; (O+A)/2 ; B ; (I+H)/2 ; P ; (K+L)/2 ; F]; 
slpsys{2,53} = [ 1  1 -2 1; -1 -1  2 6; (O+B)/2 ; C ; (J+I)/2 ; P ; (L+G)/2 ; A]; 
slpsys{2,54} = [-1  2 -1 1;  1 -2  1 6; (O+C)/2 ; D ; (K+J)/2 ; P ; (G+H)/2 ; B]; 
slpsys{2,55} = [-2  1  1 1;  2 -1 -1 6; (O+D)/2 ; E ; (L+K)/2 ; P ; (H+I)/2 ; C]; 
slpsys{2,56} = [-1 -1  2 1;  1  1 -2 6; (O+E)/2 ; F ; (G+L)/2 ; P ; (I+J)/2 ; D]; 
slpsys{2,57} = [ 1 -2  1 1; -1  2 -1 6; (O+F)/2 ; A ; (H+G)/2 ; P ; (J+K)/2 ; E]; 

%  {2112}<211-3> C2 twins:; 0.22; 1.2 Ti Zr Re]; 
slpsys{2,58} = [ 2 -1 -1 2;  2 -1 -1 -3; (J+P)/2 ; K ; F ; (F+B)/2 ; B ; I]; 
slpsys{2,59} = [ 1  1 -2 2;  1  1 -2 -3; (K+P)/2 ; L ; A ; (A+C)/2 ; C ; J]; 
slpsys{2,60} = [-1  2 -1 2; -1  2 -1 -3; (L+P)/2 ; G ; B ; (B+D)/2 ; D ; K]; 
slpsys{2,61} = [-2  1  1 2; -2  1  1 -3; (G+P)/2 ; H ; C ; (C+E)/2 ; E ; L]; 
slpsys{2,62} = [-1 -1  2 2; -1 -1  2 -3; (H+P)/2 ; I ; D ; (D+F)/2 ; F ; G]; 
slpsys{2,63} = [ 1 -2  1 2;  1 -2  1 -3; (I+P)/2 ; J ; E ; (E+A)/2 ; A ; H]; 

%   {1011}<101-2> C1 twins: 0.10; 1.1; Mg; Zr Ti]; 
slpsys{2,64} = [ 1  0 -1 1;  1  0 -1 -2; P ; L ; A ; (A+B)/2 ; B ; I ]; 
slpsys{2,65} = [ 0  1 -1 1;  0  1 -1 -2; P ; G ; B ; (B+C)/2 ; C ; J ]; 
slpsys{2,66} = [-1  1  0 1; -1  1  0 -2; P ; H ; C ; (C+D)/2 ; D ; K ]; 
slpsys{2,67} = [-1  0  1 1; -1  0  1 -2; P ; I ; D ; (D+E)/2 ; E ; L ]; 
slpsys{2,68} = [ 0 -1  1 1;  0 -1  1 -2; P ; J ; E ; (E+F)/2 ; F ; G ]; 
slpsys{2,69} = [ 1 -1  0 1;  1 -1  0 -2; P ; K ; F ; (F+A)/2 ; A ; H ]; 


% cell 3 = Cartesian slip system unit vectors, cell 4 = Schmid matrix %
for isc=1:1:nslphex 
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
     linedir(isc,:) = [isc, cross(unit_m, unit_n)];   % edge dislocation line direction unit vectors are generated and sorted
    if linedir(isc,2)<0
        linedir(isc,:) = -1.*linedir(isc,:);
    end
end
sortlinedir = sortrows(linedir,[-2,-3,-4]);
%  lots of code to set up inverse pole figure and labeling, probably has some inconsistencies
sij = [0.9581 -0.4623 -0.1893 0.698 2.1413 2.8408]/100; % for Ti from Simmons and Wang  in units of 1/GPa
Elow = 83.2;   Ehigh = 145.5;   % Elastic Ccontants of Nb from Simmons and Wang for Ti
CTE = [15.4 15.4 30.6]; % this is for Sn, not a cubic material ... such as this inconsistency

figure;  hold on; axis square; xmax = 0;
if c_a ~= 1
    trideg = 45;
    if ihex == 1
        trideg = 30;
    end
    angle = 0:1:trideg;
    xang = cosd(angle);   yang = sind(angle);
    borderx = [0 xang 0];  bordery = [0 yang 0]; xmax = 1.02; 
    axis([0 xmax 0 xmax]), % TickDir, 'out'  ???
    plot(borderx,bordery,'k-');
%     for i = 1:1:trideg+1;
%         pip = [1 (i-1)/trideg 1/c_a];   nn = norm(pip); 
%         if stereo ~=1
%             pip(3) = 0;
%         end
%         xang(i) = pip(1)/nn/(1+pip(3)/nn);
%         yang(i) = pip(2)/nn/(1+pip(3)/nn);
%     end
%     borderx = [0 xang 0];  bordery = [0 yang 0];  
%    plot(borderx,bordery,'k-');
%     for i = 1:1:trideg+1;
%         pip = [1 (i-1)/trideg c_a];   nn = norm(pip);
%         if stereo ~=1
%             pip(3) = 0;
%         end
%         xang(i) = pip(1)/nn/(1+pip(3)/nn);
%         yang(i) = pip(2)/nn/(1+pip(3)/nn);
%     end
%     borderx = [0 xang 0];  bordery = [0 yang 0];  
%    plot(borderx,bordery,'k-');
else
%     for i = 1:1:trideg+1;
%         pip = [1 (i-1)/trideg 1];   nn = norm(pip); 
%         if stereo ~=1
%             pip(3) = 0;
%         end
%         xang(i) = pip(1)/nn/(1+pip(3)/nn);
%         yang(i) = pip(2)/nn/(1+pip(3)/nn);
%     end
%     xmax = round((max(xang)*100+5)/5)/20;
%     borderx = [0 xang 0];  bordery = [0 yang 0];  
%     axis([0 xmax 0 xmax]), % TickDir, 'out'  ???
%     plot(borderx,bordery,'k-');
end

edgecolor = [0,0,1];  %  perimeter of plotting symbol
if ipfd == 2
    edgecolor = [.9,.9,0];
end
if ipfd == 1
    edgecolor = [1,0,0];
end
for p = 0:.05:1    % plot inverse pole figure gray scale symbols
    plot(.03+p*xmax*.7,.95*xmax,'o','LineWidth',4,'MarkerEdgeColor',[1 1-p p],...
                'MarkerFaceColor',[1 1 1],'MarkerSize',8)
    if c_a ~= 1
        plot(.03+p*xmax*.7,.85*xmax,'o','LineWidth',1,'MarkerEdgeColor',[p p p],...
                'MarkerFaceColor',[p p p],'MarkerSize',6)
    end
end
text(0.02,.9*xmax,'E direction, yellow (low) --> magenta (high)');
if c_a ~= 1
    text(0.02,.8*xmax,'CTE direction, black(low) --> white(high)');
end
if stereo == 0 
    text(0.02,.7*xmax,'Z projection, not stereographic');
else
    text(0.02,.7*xmax,'Stereographic projection');
end
text(0.02,.5*xmax,['IPF direction ',num2str(ipfd)]);
text(0.02,.6*xmax,['c/a ratio ',num2str(c_a)]);

%  Generate orientation matrices for each orientation in Ang
for iAng=1:1:sEul(1,1) %1;%4; %2; %  
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
    gA=gphi2*gPhi*gphi1; gAR = gA*Rotateview; % to rotate point of view in plot
    g(:,:,iAng)=gA; R(:,:,iAng)=gA';
    gsgTA = gA*ststens_n*gA'; %rotated stress tensor
    basalY(iAng,1) = abs(gA(3,2));
    
    for i = 1:1:3
        CTEm(:,i) = gA(:,i).*CTE';
    end
    % This may be a leftover from work with Sn, I don't remember what crystal structure this is for...
    % E = ( s11 * (x^4 + y^4) + s33 * z^4 + (s66 + 2*s12)*x^2*y^2 + (2*s13 + s44) * z^2 * (y^2* + x^2) )^-1*100
    for i = 1:1:3   %  This gives the CTE in the x (1), y (2), z (3) directions
%        CTExyz(iAng,i) = (CTEm(1,i)^2+CTEm(2,i)^2+CTEm(3,i)^2)^.5;
        e1 = sij(1)*(gA(1,i)^4 + gA(2,i)^4) + sij(4)*gA(3,i)^4; 
        e2 = (2*sij(2) + sij(6))*(gA(1,i)^2 * gA(2,i)^2);
        e3 = (2*sij(3) + sij(5))*gA(3,i)^2 * (gA(1,i)^2 + gA(2,i)^2);
        EA(iAng,i)=1./(e1+e2+e3);   
    end
%    CTEgray = (CTExyz(iAng,CTEd)-CTE(1)*0.999999)/(CTE(3)-CTE(1)*0.999999);%+1e-14;
    Eb = (EA(iAng,iEd)-Elow)/(Ehigh-Elow);
%     XPgray = (Ang(iAng,1)-XPgl)/(XPgh-XPgl);
%     if ipfd == 2 
%         XPgray = 1-XPgray;
%     end
%    gray = [CTEgray CTEgray CTEgray];   %XPgray XPgray XPgray
    pv = gA(:,ipfd)';    pv = abs(pv);
    if c_a == 1
        x(1) = median(pv);
        x(2) = min(pv);
        x(3) = max(pv);
    else
        x(1) = max(pv(1),pv(2));
        x(2) = min(pv(1),pv(2));
        x(3) = pv(3);
    end
    if ihex == 1
        ang = atand(x(2)/x(1));
        % [iAng x ang]
        if ang > 30 % tangent of 30 deg
            nang = 30-(ang - 30);
            radius = (x(1)^2+x(2)^2)^.5;
            x(1) = cosd(nang)*radius;
            x(2) = sind(nang)*radius;
        end            
    end
    if stereo == 1
        plot(x(1)/(1+x(3)), x(2)/(1+x(3)),'o','LineWidth',2,'MarkerEdgeColor',...
                edgecolor,'MarkerFaceColor',[1 1-Eb Eb],'MarkerSize',14)
%         plot(x(1)/(1+x(3)), x(2)/(1+x(3)),'o','LineWidth',1,'MarkerEdgeColor',...
%                  [1 1 1],'MarkerFaceColor',gray,'MarkerSize',6)
        text(x(1)/(1+x(3)), x(2)/(1+x(3))+.037,num2str(iAng))
    else
%         plot(x(1),x(2),'o','LineWidth',1,'MarkerEdgeColor',edgecolor,...
%                 'MarkerFaceColor',gray,'MarkerSize',14)
        plot(x(1),x(2),'o','LineWidth',1,'MarkerEdgeColor',edgecolor,...
                 'MarkerFaceColor',[1 1-Eb Eb],'MarkerSize',8)
        text(x(1),x(2)+.037,num2str(iAng))
    end
    
%  Compute Schmid Factors, followed by grain slip system; slpsys{4 = Schmid matrix}            
    for isc=1:1:nslphex 
        computeSFA(isc,1,iAng)=isc;   %This variable will have slip system number, Schmid Factor, plane and Burgers in sample coord syst, and hkl,uvw
        computeSFA(isc,2,iAng)=0.;
        for i=1:1:3
            for j=1:1:3
                computeSFA(isc,2,iAng)=computeSFA(isc,2,iAng)+gsgTA(i,j)*slpsys{4,isc}(i,j);
            end
        end               % Abs(Schmid factor) is in column after signed Schmid factor
        if isc>= 46 && computeSFA(isc,2,iAng)<0 
            computeSFA(isc,2,iAng) = 0.001*computeSFA(isc,2,iAng) ;    % this is to prevent anti-twin shears from being seriously considered later
        end
        rot_nA = slpsys{3,isc}(1,:)*gAR;
        rot_bA = slpsys{3,isc}(2,:)*gAR;
        computeSFA(isc,3:17,iAng)=[abs(computeSFA(isc,2,iAng)),rot_nA, rot_bA, slpsys{2,isc}(1,:),slpsys{2,isc}(2,:)];
        rot_p1 = slpsys{3,isc}(3,:)*gAR;
        rot_p2 = slpsys{3,isc}(4,:)*gAR;
        rot_p3 = slpsys{3,isc}(5,:)*gAR;
        rot_p4 = slpsys{3,isc}(6,:)*gAR;
        rot_p5 = slpsys{3,isc}(7,:)*gAR;
        rot_p6 = slpsys{3,isc}(8,:)*gAR;
% Schmid factors (1-3), rotated plane normal (4-6), Computed rotated Burgers vector (7-9),  
% , plane trace on Z surface (10-12) Computed rotated position vectors to points p1-p4 (13-24)
        Schm_labvecA(isc,:) = [computeSFA(isc,1:3,iAng), rot_nA, rot_bA, cross(rot_nA',[0,0,1]), rot_p1, rot_p2, rot_p3, rot_p4, rot_p5, rot_p6];
    end                             %    1-3         4-6      7-9               10-12        13-15   16-18   19-21   22-24   25-27   28-30

% sortmv contains list from high to low Schmid factor, with associated infomation to draw slip system in unit cell
% useful plotting unit cell vectors will sort to bottom row
    Schm_labvecA(nslphex+1,:) = [0 1 -1 [1 0 0]*gAR [0 1 0]*gAR [0 0 1]*gAR 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
    sortmv(:,:,iAng) = sortrows(Schm_labvecA,-3);
end
%%  plot the image of unit cell, slip vectors, planes, plane normals, and plane traces
%  Then choose what row (orientation) to analyze, 
%  and set iAng to this orientation (row), and run this second cell.    X down and Y to right!!!
iAng = 1;   
ptpl = 1;  % 0 = don't plot plane traces,
pta123o = 1;  % 1 = a1 a2 a3 coordinates and origin (open circle), 
ptss = 3; % 1 slip planes, 2 planes and directions, 3 planes, directions, normals,  4 directions only, 
nplts = nslphex; %48; %  number of slip sorted systems to plot
%  no further user input below here

%  Dashed lines give plane traces, (colors in groups based on slip
%  families) Shorter plane traces imply that the slip plane is nearly
%  parallel to page, longer traces imply that the plane is highly inclined.
%  dotted red, green, blue lines give the x,y,z edges of unit cell.  
%  Turquiose/Teal line shows Burgers vector direction, with direction
%  away from the ball end (--!--> adjusted for the sign of the Schmid factor <--!--).
%  Slip plane is shaded light when the plane normal has out-of-page
%  component, darker when into the page.  

%  Strategy:  First extract useful vectors from slip system information to draw the hexagonal prisms
%  positions in sortmv  p1:13-15  p2:16-18  p3:19-21  p4:22-24  p5:25-27  p6:28-30 
%  positions in pln    p1:4-6   p2:7-9   p3:10-12  p4:13-15  p5:16-18  p6:18-21
for isc = 1:1:nslphex   
    if sortmv(isc,1,iAng) == 1               % locate the two basal planes
        pln(1,4:21) = sortmv(isc,13:30,iAng);  % bottom basal plane
        pln(2,4:21) = sortmv(isc,13:30,iAng);  % top basal plane
        rotc = sortmv(isc,4:6,iAng)*c_a;       % basal plane normal * c/a
        for j = 4:3:19
            pln(2,j:j+2) = pln(1,j:j+2) + rotc;    % move top plane up by a unit of c
        end
        a1 = sortmv(isc,7:9,iAng);   %  locate a1 using SS1
    elseif sortmv(isc,1,iAng) == 2
        a2 = sortmv(isc,7:9,iAng);   %  locate a2 using SS2
    elseif sortmv(isc,1,iAng) == 3
        a3 = sortmv(isc,7:9,iAng);   %  locate a3 using SS3
    end
end
for isc = 1:1:nslphex
    if sortmv(isc,1,iAng) == 4        %  locate two prism planes on opposite sides using SS4
    	pln(3,4:21) = sortmv(isc,13:30,iAng);
        for j = 13:3:28
            pln(4,j-9:j-7) = sortmv(isc,j:j+2,iAng) + a2 - a3;
        end
    elseif sortmv(isc,1,iAng) == 5    %  locate two prism planes on opposite sides using SS5
    	pln(5,4:21) = sortmv(isc,13:30,iAng);
     	for j = 13:3:28
            pln(6,j-9:j-7) = sortmv(isc,j:j+2,iAng) + a3 - a1;
        end
    elseif sortmv(isc,1,iAng) == 6    %  locate two prism planes on opposite sides using SS6
    	pln(7,4:21) = sortmv(isc,13:30,iAng);
     	for j = 13:3:28
            pln(8,j-9:j-7) = sortmv(isc,j:j+2,iAng) + a1 - a2;
        end
    end
end
for j = 1:1:2    % Find z elevation of basal planes
    for k = 1:1:3
        pln(j,k) = (pln(j,3+k)+pln(j,6+k)+pln(j,9+k)+pln(j,12+k)+pln(j,15+k)+pln(j,18+k))/6;
    end
end
center = (pln(1,1:3)+pln(2,1:3))/2; 
for j = 3:1:8  % Find z elevation of prism planes
    pln(j,3) = (pln(j,12)+pln(j,15)+pln(j,18)+pln(j,21))/4;
end
sortpln = sortrows(pln,-3);
minx = 0; miny = 0; minz = 0; maxx = 0; maxy = 0; maxz = 0;
for j = 1:1:8 % assemble vectors for plotting faces of hex prism
    prsxplt(j,1:7) = [sortpln(j,4) sortpln(j,7) sortpln(j,10) sortpln(j,13) sortpln(j,16) sortpln(j,19) sortpln(j,4)];
    minx = min(minx,min(prsxplt(j,:))); maxx = max(maxx,max(prsxplt(j,:)));
    prsyplt(j,1:7) = [sortpln(j,5) sortpln(j,8) sortpln(j,11) sortpln(j,14) sortpln(j,17) sortpln(j,20) sortpln(j,5)];
    miny = min(miny,min(prsyplt(j,:))); maxy = max(maxy,max(prsyplt(j,:)));
    prszplt(j,1:7) = [sortpln(j,6) sortpln(j,9) sortpln(j,12) sortpln(j,15) sortpln(j,18) sortpln(j,21) sortpln(j,6)];
    minz = min(minz,min(prszplt(j,:))); maxz = max(maxz,max(prszplt(j,:)));
end

ipl = -8; %  Strategy: Next, start isc loop for plotting slip systems
for isc = 1:1:nplts  %nslphex % change  to smaller number to plot fewer nslphex 
    if ipl-isc==-9 % eight plots on a page
        figure
        ipl=ipl+8;
    end
    subplot(2,4,isc-ipl)
    hold on
    sp1 = sortmv(isc,13:15,iAng);       % beginning of Burgers vector
    sp2 = sortmv(isc,16:18,iAng);       % extract plotted points on perimeter of the slip plane
    sp3 = sortmv(isc,19:21,iAng);
    sp4 = sortmv(isc,22:24,iAng);
    sp5 = sortmv(isc,25:27,iAng);
    sp6 = sortmv(isc,28:30,iAng);
    spx = [sp1(1) sp2(1) sp3(1) sp4(1) sp5(1) sp6(1) sp1(1)];
    spy = [sp1(2) sp2(2) sp3(2) sp4(2) sp5(2) sp6(2) sp1(2)];
    ssn = sortmv(isc,1,iAng);          % slip system number
    Sf = sortmv(isc,2,iAng);           % Schmid factor
    Sfs = 1;
    if Sf < 0
        Sfs = -1;
    end
    n = [0 0 0 sortmv(isc,4:6,iAng)];  % plane normal
    b = [sp1 sp4]; % p1+sortmv(isc,7:9,iAng)];  % Burgers vector
    nvec = sortmv(isc,4:6,iAng);
    bvec = sortmv(isc,7:9,iAng);
    pt = sortmv(isc,10:12,iAng);       % plane trace
    minx = min(minx, sp1(1)+n(4));
    maxx = max(maxx, sp1(1)+n(4));   % find appropriate range of x and y for plot
    miny = min(miny, sp1(2)+n(5));
    maxy = max(maxy, sp1(2)+n(5));
    midx = (minx+maxx)/2;
    midy = (miny+maxy)/2;
    del = 2;

% These plots will match TSL with X down !!!!   Plotting starts...
    axis square
    set(gca ,'ycolor' ,'w'); set(gca ,'xcolor' ,'w');  % make axes white for ease in later arranging.
    axis([midy-del midy+del -midx-del -midx+del])
    if pta123o > 0
        if Ang(iAng,4) < 90 % make the 3 coordinate axes visible below slip planes
            plot([0 a1(2)],-[0 a1(1)], ':', 'Linewidth',3,'Color',[1 0 .2]);% plot x = red
            plot([0 a2(2)],-[0 a2(1)], ':', 'Linewidth',3,'Color',[.6 .8 0]);% plot y = green-gold
            plot([0 a3(2)],-[0 a3(1)], ':', 'Linewidth',3,'Color',[0 0 1]);% plot z = blue
        end
    end
    if ipl-isc==-1
       %text(min_x, min_y-1, 'Eulers = ',num2str(Ang(iAng,3)), num2str(Ang(iAng,4)), num2str(Ang(iAng,5)));
    end
    if ptss <= 3
        if sortmv(isc,6,iAng) > 0   % is k component of slip plane normal positive or negative?
            fill(spy,-spx, [.8 .8 .65])  % slip plane filled warm gray
    %        plot([n(2) n(5)], -[n(1) n(4)],'Linewidth',3,'Color',[.8 .8 .65]);
        else                % slip plane filled cool gray if normal has neg z component
            fill(spy,-spx, [.65 .65 .7])  
    %        plot([n(2) n(5)], -[n(1) n(4)],'Linewidth',3,'Color',[.65 .65 .7]);
        end
        if n(6) > 0
            pncolor = [0 0 0];
        else
            pncolor = [.5 .5 .5];
        end
        if sortmv(isc,6,iAng) > 0
            Bvcolor = [0 .7 .7];
            if ssn >= 58
                Bvcolor = [.1 .6 0];
            end
            if ssn >= 46 && ssn < 58
                Bvcolor = [1 .6 0];
            end
        else
            Bvcolor = [0 1 1];
            if ssn >= 58
                Bvcolor = [.3 .9 0];
            end
            if ssn >= 46 && ssn < 58
                Bvcolor = [1 .8 0];
            end
        end
        Sfs = 1;
        if ssn < 46   % will reverse the sign of Schmid factor for dislocations
        end
        if ptss > 1
            if Sf > 0    % plot Burgers vector direction
                if ssn >= 46    % this is for twins - the Burgers vector length is shown to be 1/2 of the usual length in the unit cell 
                    plot(b(2),-b(1),'.','MarkerSize', 24, 'Color', Bvcolor)
                    plot([b(2) (b(2)+b(5))/2],-([b(1) (b(1)+b(4))/2]),'Linewidth',4,'Color',Bvcolor)
                else
                    plot(b(2),-b(1),'.','MarkerSize', 24, 'Color', Bvcolor)
                    plot([b(2) b(5)],-[b(1) b(4)],'Linewidth',4,'Color',Bvcolor)
    %                 quiver(p1(1),p1(2),dp(1),dp(2),0,'Linewidth',2,'Color',Bvcolor)
                end
            else         % plot Burgers vector in opposite direction
                if ssn >= 46    % this is for twins - the Burgers vector length is shown to be 1/2 of the usual length in the unit cell 
                    plot(b(2),-b(1),'.','MarkerSize', 24, 'Color', Bvcolor)
                    plot([b(2) (2*b(2)+b(5))/3],-([b(1) (2*b(1)+b(4))/3]),'Linewidth',4,'Color',Bvcolor)
                else
                    plot(b(5),-b(4),'.','MarkerSize', 24, 'Color', Bvcolor)
                    plot([b(5) b(2)],-[b(4) b(1)],'Linewidth',4,'Color',Bvcolor)
    %                quiver(p2(1),p2(2),-dp(1),-dp(2),0,'Linewidth',4,'Color',Bvcolor)
                end
            end
            if ptss ==3
                plot([b(2) (b(2)+n(5))],-([b(1) (b(1)+n(4))]),'Linewidth',4,'Color',pncolor)
            end
        end
    end
    for j = 1:1:4 % plot the 4 top most surface prisms of the hex cell that have the highest z elevation
        plot(prsyplt(j,:),-prsxplt(j,:), 'Linewidth',2,'Color',[.0 .0 .0]);
    end
    if pta123o > 0
        if Ang(iAng,4) > 90 % make the 3 coordinate axes visible above slip planes
            plot([0 a1(2)],-[0 a1(1)], ':', 'Linewidth',3,'Color',[1 0 .3]);% plot x = red
            plot([0 a2(2)],-[0 a2(1)], ':', 'Linewidth',3,'Color',[.5 .6 0]);% plot y = green-gold
            plot([0 a3(2)],-[0 a3(1)], ':', 'Linewidth',3,'Color',[0 0 1]);% plot z = blue
        end
    end
    if ptpl == 1   % plot plane traces
        if ssn>=58  % compression twin plane traces   green
            plot([-pt(2) pt(2)],-[-pt(1) pt(1)],'--','Linewidth',3,'Color',[.2 .8 0]) 
        elseif ssn>45 && ssn<58  % extension twin plane traces   orange
            plot([-pt(2) pt(2)],-[-pt(1) pt(1)],'--','Linewidth',3,'Color',[1 .6 0]) 
        elseif ssn>15 && ssn<46   % <c+a> plane traces   green-gold
            plot([-pt(2) pt(2)],-[-pt(1) pt(1)],'--','Linewidth',3,'Color',[.95 .85 0]) 
        elseif ssn<16 && ssn>9   % pyr <a>  green
            plot([-pt(2) pt(2)],-[-pt(1) pt(1)],'--','Linewidth',3,'Color',[0 .9 .5]) 
        elseif ssn<10 && ssn>3   % prism <a>  red
            plot([-pt(2) pt(2)],-[-pt(1) pt(1)],'--','Linewidth',3,'Color',[1 .2 0]) 
        else          %  {medium slip systems}  blue 
            plot([-pt(2) pt(2)],-[-pt(1) pt(1)],'--','Linewidth',3,'Color',[0 0 1]) 
        end                %---->  NOTE that Schmid factor vector is plotted in correct direction, 
    end
    if pta123o == 1
        plot(0,0,'ko'); 
    end
%---->  Burgers vector is labeled and plotted with consistently signed b vector direction.
    if ptss >= 1
        title({['ssn' num2str(ssn) ' n' mat2str(slpsys{2,ssn}(1,:))  mat2str(Sfs*slpsys{2,ssn}(2,:)) 'b'],... 
        ['Or-' num2str(iAng) ' m' num2str(isc) ' = ' num2str(Sfs*Sf, 4) ],...
        ['Eulers = ', mat2str(Ang(iAng,3:5),3)]}) 
    end
end
%          ['n = ', mat2str(nvec,3)],...
%          ['b = ', mat2str(bvec,3)],...



end