%% Function to evaluate Slip Transfer Parameter (Luster-Morris) given Euler angles and 
%% Stress Tensors in sample coordinate system
% This is based on Dr. Bieler's code for Plotting mprime,fip and Schmid
% factor. The difference here is that this is customised to read data from
% APS Beamline 1
% Arranges grains in a 2x2 layout
function MPrimeFun(EulerTL,EulerTR,EulerBL,EulerBR,STens,ScTol,GBNumChoice)
%Inputs: EulerBL:Euler angle of bottom left grain, EulerTL: Euler angle of
%top left grain, EulerTR: Euler angle of top right grain
% Stens:Stress Tensor(1x6) tuple
%SchTol:Schmid Factor tolerance [Low,High]
% First Step: Generate Grain ID and RCB Arrays
[IDgr,RBdy]=GIDRCB_Gen(EulerTL,EulerTR,EulerBL,EulerBR);

dRBdy = size(RBdy);  % This is a reconstructed grain boundary file 
dIDgr = size(IDgr);  % This is a default type 2 grain file, both are needed.
nslphex = 39;    c_a_hex = 1.59;  % 1.587 for pure Ti ... not turning on compression twinning in Hexagonal
nslpbcc = 24;    c_a_bcc = 1.0;    % could differ for metastable phases... not turning on 123 slip 
nslpfcc = 12;    c_a_fcc = 1.0;    % ~1.02 for TiAl%  
nslpbct = 32;    c_a_bct = 0.5456; %  for Sn   
nslp = [nslphex, nslpbcc, nslpfcc, nslpbct]; % number of slip systems used for phases 1, 2, 3, 4
c_a = [c_a_hex, c_a_bcc, c_a_fcc, c_a_bct];

% If your data is SINGLE PHASE, then you must put the correct phase number into the variable one_ss, HERE
one_ss = 1; %  e.g. for BCC, the if statement below will set phase  = 2 (3 for FCC) in IDgr file column 10:
numPhases = max(IDgr(:,10));   % check the number of phases in the dataset 
if numPhases == 0            % and correct, if needed 
    IDgr(:,10) = one_ss;     % single phase, ---!!! set to 1 for HEX or 2 for BCC above !!!---
    for i = 1:1:4             % else if it's two phase, do nothing, if phase 1 is hex, 2 is bcc.
        if i ~= one_ss
            nslp(i) = 0;
        end
    end
end                

schmid_tolH = ScTol(2);%0.25; % higher schmid tolerance value
schmid_tolL = ScTol(1);%0.0; % lower schmid tolerance value 
hkl = 1;           % flag used to decide whether to adjust first euler angle for various reasons...  see below

% Stress tensor is defined using TSL convensions with x down !!!  put the one you want last
% sigma = [1,0,0; 0,0,0; 0,0,0];  sigma = [0,1,0; 1,0,0; 0,0,0];  sigma = [0,0,0; 0,1,0; 0,0,0]; 
sigma = [STens(1),STens(4),STens(5); 
         STens(4),STens(2),STens(6); 
         STens(5),STens(6),STens(3)]; 

nsten = 1;
%ststens = [0 0 0 ; 0 1 0 ; 0 0 0];  tension in Y  %  [.9 .01 .02 ; .01 0 0 ; .02 0 0]
for i = 1:1:nsten
    str2 = sigma(:,:,i)*sigma(:,:,i)';
    ststens_mag = (str2(1,1)+str2(2,2)+str2(3,3))^.5;
    sigma_n(:,:,i) = sigma(:,:,i)/ststens_mag;   % normalized stress tensor to get generalized Schmid factor
    sigma_v(:,i) = [sigma(1,1) sigma(2,2) sigma(3,3)]';  % vectorized version of trace
end
% E(X,Y,Z) will be calculated from [S11 S12 S13 S33 S44 S66] in 1/GPa; 
% hexagonal stiffness chosen:
sij(1,:) = [0.9581 -0.4623 -0.1893 0.698 2.1413 2.408]/100; % for Ti from Simmons and Wang  in units of 1/GPa 
% cubic stiffness chosen:
sij(2,:) = [0.6862 -0.2581 -0.2581 0.6862 1.2123 1.2123]/100; % for Ta from Simmons and Wang  in units of 1/GPa 
%  disp('E(X,Y,Z) in GPa; S11 S12 S13 S33 S44 S66  210 Rayne, J.A. and B.S. Chandrasekhar, 
%  Elastic Ccontants of  beta tin from 4.2K to 300K, Phys Rev. 118, 1545-49, 1960
sij(4,:) = [4.3627 -3.3893 -0.394 1.4501 4.5393 4.1667]/100; % for Sn in units of 1/GPa 
% Ti-6Al, Ti-15Cr, alpha/beta in Ti6242 from J. Kim and S.I. Rokhlin, J. Acoust. Soc. Am. 126-6 dec 2009
% Sij [2.56 -1.14 -1.14 2.56 2.26 2.26]/100; % for beta Ti-15Cr 
% Sij [3.12 -1.42 -1.42 3.12 1.82 1.82]/100; % for beta in Ti-6242 
% Sij [1.182 -0.5424 -0.2686 0.839 2.463 3.448]/100; % for Ti-6Al 
% Sij [1.063 -0.4972 -0.2010 0.756 2.053 3.120]/100; % for alpha in Ti-6242 
% Sij [0.9581 -0.4623 -0.1893 0.698 2.1413 2.408]/100; % for Ti from Simmons and Wang 
% Sij [0.6874 -0.2481 -0.2481 0.6874 3.413 3.413]/100; % for Nb from Simmons and Wang  
% Sij [0.6862 -0.2581 -0.2581 0.6862 1.2123 1.2123]/100; % for Ta from Simmons and Wang  
% 


% %  Set up vectors useful for plotting unit cells with slip systems


O = [ 0  0  0 0];   %          (I)         (H)
A = [ 2 -1 -1 0]/3; %            C ------- B 
B = [ 1  1 -2 0]/3; %          /  \       / \
C = [-1  2 -1 0]/3; %         /    a2   /    \
D = [-2  1  1 0]/3; %        /       \ /      \
E = [-1 -1  2 0]/3; %    (J)D ------ O(P)-a1-> A(G) --> x
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
sshex(:,:,1) = [0 0 0 1;  2 -1 -1 0; D ; E ; F ; A ; B ; C ]; 
sshex(:,:,2) = [0 0 0 1; -1  2 -1 0; F ; A ; B ; C ; D ; E ]; 
sshex(:,:,3) = [0 0 0 1; -1 -1  2 0; B ; C ; D ; E ; F ; A ]; 
ibas = 1;
fbas = 3;
%  prism <a>-glide:Ti Zr RE; Be Re Mg
sshex(:,:,4) = [ 0  1 -1 0;  2 -1 -1 0; E ; E ; E ; F ; L ; K ];
sshex(:,:,5) = [-1  0  1 0; -1  2 -1 0; A ; A ; A ; B ; H ; G ];
sshex(:,:,6) = [ 1 -1  0 0; -1 -1  2 0; C ; C ; C ; D ; J ; I ];
iprs = 4;
fprs = 6;
% prism <aa>
sshex(:,:,7) = [ 2 -1 -1 0;  0  1 -1 0; F; F; F; B ; H ; L ];
sshex(:,:,8) = [-1  2 -1 0; -1  0  1 0; B; B; B; D ; J ; H ];
sshex(:,:,9) = [-1 -1  2 0;  1 -1  0 0; D; D; D; F ; L ; J ];
i2prs = 7;
f2prs = 9;
%  pyramidal <a>-glide --**-- CORRECTED --**--
sshex(:,:,10) = [ 0  1 -1 1;  2 -1 -1 0; C ; C ; C ; B ; G ; J ];
sshex(:,:,11) = [-1  0  1 1; -1  2 -1 0; E ; E ; E ; D ; I ; L ];
sshex(:,:,12) = [ 1 -1  0 1; -1 -1  2 0; A ; A ; A ; F ; K ; H ];
sshex(:,:,13) = [-1  1  0 1;  1  1 -2 0; D ; D ; D ; C ; H ; K ];
sshex(:,:,14) = [ 0 -1  1 1; -2  1  1 0; F ; F ; F ; E ; J ; G ];
sshex(:,:,15) = [ 1  0 -1 1;  1 -2  1 0; B ; B ; B ; A ; L ; I ];
ipyra = 10;
fpyra = 15;
%  pyramidal <c+a>-glide:; all?
sshex(:,:,16) = [-1  1  0 1;  2 -1 -1 3; D ; D ; K ; P ; H ; C ]; 
sshex(:,:,17) = [-1  1  0 1;  1 -2  1 3; C ; D ; K ; P ; H ; C ];
sshex(:,:,18) = [ 1  0 -1 1; -1 -1  2 3; B ; B ; I ; P ; L ; A ];
sshex(:,:,19) = [ 1  0 -1 1; -2  1  1 3; A ; B ; I ; P ; L ; A ];
sshex(:,:,20) = [ 0 -1  1 1; -1  2 -1 3; F ; F ; G ; P ; J ; E ];
sshex(:,:,21) = [ 0 -1  1 1;  1  1 -2 3; E ; F ; G ; P ; J ; E ];
sshex(:,:,22) = [ 1 -1  0 1; -2  1  1 3; A ; A ; H ; P ; K ; F ];
sshex(:,:,23) = [ 1 -1  0 1; -1  2 -1 3; F ; A ; H ; P ; K ; F ];
sshex(:,:,24) = [-1  0  1 1;  1  1 -2 3; E ; E ; L ; P ; I ; D ];
sshex(:,:,25) = [-1  0  1 1;  2 -1 -1 3; D ; E ; L ; P ; I ; D ]; 
sshex(:,:,26) = [ 0  1 -1 1;  1 -2  1 3; C ; C ; J ; P ; G ; B ];
sshex(:,:,27) = [ 0  1 -1 1; -1 -1  2 3; B ; C ; J ; P ; G ; B ];
ipyrc = 16;
fpyrc = 27;
%  pyramidal <c+a>-2nd order glide
sshex(:,:,28) = [-2  1  1 2;  2 -1 -1 3; (O+D)/2 ; E ; L ; (P+G)/2 ; H ; C]; 
sshex(:,:,29) = [ 1 -2  1 2; -1  2 -1 3; (O+F)/2 ; A ; H ; (P+I)/2 ; J ; E]; 
sshex(:,:,30) = [ 1  1 -2 2; -1 -1  2 3; (O+B)/2 ; C ; J ; (P+K)/2 ; L ; A]; 
sshex(:,:,31) = [ 2 -1 -1 2; -2  1  1 3; (O+A)/2 ; B ; I ; (P+J)/2 ; K ; F]; 
sshex(:,:,32) = [-1  2 -1 2;  1 -2  1 3; (O+C)/2 ; D ; K ; (P+L)/2 ; G ; B]; 
sshex(:,:,33) = [-1 -1  2 2;  1  1 -2 3; (O+E)/2 ; F ; G ; (P+H)/2 ; I ; D]; 
i2pyrc = 28;
f2pyrc = 33;
% *** Twin directions are opposite in Christian and Mahajan, and are not correcte to to be consistent with them
% FROM Kocks SXHEX   plane   direction   1st and 4th point is Burgers vector, order of C1 differs from Kock's file
%  {1012}<1011> T1 twins 0.17; -1.3  twins: all  Twin Vector must go in the
%  sense of shear, opposite C&M sense.
sshex(:,:,34) = [-1  1  0 2;  1 -1  0 1; C ; D ; L ; G ; G ; G ]; 
sshex(:,:,35) = [ 1  0 -1 2; -1  0  1 1; A ; B ; J ; K ; K ; K ]; 
sshex(:,:,36) = [ 0 -1  1 2;  0  1 -1 1; E ; F ; H ; I ; I ; I ]; 
sshex(:,:,37) = [ 1 -1  0 2; -1  1  0 1; F ; A ; I ; J ; J ; J ]; 
sshex(:,:,38) = [-1  0  1 2;  1  0 -1 1; D ; E ; G ; H ; H ; H ]; 
sshex(:,:,39) = [ 0  1 -1 2;  0 -1  1 1; B ; C ; K ; L ; L ; L ]; 
iT1 = 34;
fT1 = 39;
%  {2111}<2116> T2 twins: 0.63;  -0.4; Ti Zr Re RE]; Also does not follow C&M definition for shear direction
sshex(:,:,40) = [-2  1  1 1;  2 -1 -1 6; (O+D)/2 ; E ; (L+K)/2 ; P ; (H+I)/2 ; C]; 
sshex(:,:,41) = [ 1 -2  1 1; -1  2 -1 6; (O+F)/2 ; A ; (H+G)/2 ; P ; (J+K)/2 ; E]; 
sshex(:,:,42) = [ 1  1 -2 1; -1 -1  2 6; (O+B)/2 ; C ; (J+I)/2 ; P ; (L+G)/2 ; A]; 
sshex(:,:,43) = [ 2 -1 -1 1; -2  1  1 6; (O+A)/2 ; B ; (I+H)/2 ; P ; (K+L)/2 ; F]; 
sshex(:,:,44) = [-1  2 -1 1;  1 -2  1 6; (O+C)/2 ; D ; (K+J)/2 ; P ; (G+H)/2 ; B]; 
sshex(:,:,45) = [-1 -1  2 1;  1  1 -2 6; (O+E)/2 ; F ; (G+L)/2 ; P ; (I+J)/2 ; D]; 
iT2 = 40;
fT2 = 45;
%   {1011}<101-2> C1 twins: 0.10; 1.1; Mg; Zr Ti]; agrees with C&M
sshex(:,:,46) = [-1  1  0 1; -1  1  0 -2; P ; H ; C ; (C+D)/2 ; D ; K ]; 
sshex(:,:,47) = [ 1  0 -1 1;  1  0 -1 -2; P ; L ; A ; (A+B)/2 ; B ; I ]; 
sshex(:,:,48) = [ 0 -1  1 1;  0 -1  1 -2; P ; J ; E ; (E+F)/2 ; F ; G ]; 
sshex(:,:,49) = [ 1 -1  0 1;  1 -1  0 -2; P ; K ; F ; (F+A)/2 ; A ; H ]; 
sshex(:,:,50) = [-1  0  1 1; -1  0  1 -2; P ; I ; D ; (D+E)/2 ; E ; L ]; 
sshex(:,:,51) = [ 0  1 -1 1;  0  1 -1 -2; P ; G ; B ; (B+C)/2 ; C ; J ]; 
iC1 = 46;
fC1 = 51;
%  {2112}<211-3> C2 twins:; 0.22; 1.2 Ti Zr Re]; agrees with C&M
sshex(:,:,52) = [ 2 -1 -1 2;  2 -1 -1 -3; (J+P)/2 ; K ; F ; (F+B)/2 ; B ; I]; 
sshex(:,:,53) = [-1  2 -1 2; -1  2 -1 -3; (L+P)/2 ; G ; B ; (B+D)/2 ; D ; K]; 
sshex(:,:,54) = [-1 -1  2 2; -1 -1  2 -3; (H+P)/2 ; I ; D ; (D+F)/2 ; F ; G]; 
sshex(:,:,56) = [-2  1  1 2; -2  1  1 -3; (G+P)/2 ; H ; C ; (C+E)/2 ; E ; L]; 
sshex(:,:,57) = [ 1 -2  1 2;  1 -2  1 -3; (I+P)/2 ; J ; E ; (E+A)/2 ; A ; H]; 
sshex(:,:,58) = [ 1  1 -2 2;  1  1 -2 -3; (K+P)/2 ; L ; A ; (A+C)/2 ; C ; J]; 
iC2 = 52;
fC2 = 58;

%   BCC            plane   direction   -------------------------------------
ssbcc(:,:,1) =  [ 1 -1  0 ;  1  1  1 ; 0 0 0 ; 1 1 0 ; 1 1 1 ; 0 0 1 ; 0 0 1 ; 0 0 1]; %ok
ssbcc(:,:,2) =  [-1  0  1 ;  1  1  1 ; 0 0 0 ; 1 0 1 ; 1 1 1 ; 0 1 0 ; 0 1 0 ; 0 1 0]; %ok
ssbcc(:,:,3) =  [ 0 -1  1 ;  1  1  1 ; 0 0 0 ; 1 0 0 ; 1 1 1 ; 0 1 1 ; 0 1 1 ; 0 1 1]; %ok

ssbcc(:,:,4) =  [ 1  1  0 ; -1  1  1 ; 1 0 0 ; 0 1 0 ; 0 1 1 ; 1 0 1 ; 1 0 1 ; 1 0 1]; %ok
ssbcc(:,:,5) =  [ 1  0  1 ; -1  1  1 ; 1 0 0 ; 1 1 0 ; 0 1 1 ; 0 0 1 ; 0 0 1 ; 0 0 1]; %ok
ssbcc(:,:,6) =  [ 0 -1  1 ; -1  1  1 ; 1 0 0 ; 1 1 1 ; 0 1 1 ; 0 0 0 ; 0 0 0 ; 0 0 0]; %ok

ssbcc(:,:,7) =  [ 1 -1  0 ; -1 -1  1 ; 1 1 0 ; 0 0 0 ; 0 0 1 ; 1 1 1 ; 1 1 1 ; 1 1 1]; %ok
ssbcc(:,:,8) =  [ 1  0  1 ; -1 -1  1 ; 0 0 1 ; 1 0 0 ; 1 1 0 ; 0 1 1 ; 0 1 1 ; 0 1 1]; %ok
ssbcc(:,:,9) =  [ 0  1  1 ; -1 -1  1 ; 1 1 0 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 1 0 1 ; 1 0 1]; %ok

ssbcc(:,:,10) = [ 1  1  0 ;  1 -1  1 ; 0 1 0 ; 0 1 1 ; 1 0 1 ; 1 0 0 ; 1 0 0 ; 1 0 0]; %ok
ssbcc(:,:,11) = [-1  0  1 ;  1 -1  1 ; 0 1 0 ; 0 0 0 ; 1 0 1 ; 1 1 1 ; 1 1 1 ; 1 1 1]; %ok
ssbcc(:,:,12) = [ 0  1  1 ;  1 -1  1 ; 0 1 0 ; 0 0 1 ; 1 0 1 ; 1 1 0 ; 1 1 0 ; 1 1 0]; %ok
i110 = 1;
f110 = 12;

% Mode 2,   plane direction, Define four points in the plane 
ssbcc(:,:,13) = [-1 -1  2 ;  1  1  1 ; 0 0 0 ; 1 0 0.5 ; 1 1 1 ; 0 1 0.5 ; 0 1 0.5 ; 0 1 0.5]; %ok
ssbcc(:,:,14) = [ 1 -2  1 ;  1  1  1 ; 0 0 0 ; 1 0.5 0 ; 1 1 1 ; 0 0.5 1 ; 0 0.5 1 ; 0 0.5 1]; %ok
ssbcc(:,:,15) = [-2  1  1 ;  1  1  1 ; 0 0 0 ; 0.5 0 1 ; 1 1 1 ; 0.5 1 0 ; 0.5 1 0 ; 0.5 1 0]; %ok

ssbcc(:,:,16) = [ 1 -1  2 ; -1  1  1 ; 1 0 0 ; 1 1 0.5 ; 0 1 1 ; 0 0 0.5 ; 0 0 0.5 ; 0 0 0.5]; %ok
ssbcc(:,:,17) = [-1 -2  1 ; -1  1  1 ; 1 0 0 ; 1 0.5 1 ; 0 1 1 ; 0 0.5 0 ; 0 0.5 0 ; 0 0.5 0]; %ok
ssbcc(:,:,18) = [ 2  1  1 ; -1  1  1 ; 1 0 0 ; 0.5 1 0 ; 0 1 1 ; 0.5 0 1 ; 0.5 0 1 ; 0.5 0 1]; %ok  

ssbcc(:,:,19) = [ 1  1  2 ; -1 -1  1 ; 1 1 0 ; 0 1 0.5 ; 0 0 1 ; 1 0 0.5 ; 1 0 0.5 ; 1 0 0.5]; %ok
ssbcc(:,:,20) = [-1  2  1 ; -1 -1  1 ; 1 1 0 ; 0 0.5 0 ; 0 0 1 ; 1 0.5 1 ; 1 0.5 1 ; 1 0.5 1]; %ok
ssbcc(:,:,21) = [ 2 -1  1 ; -1 -1  1 ; 1 1 0 ; 0.5 1 1 ; 0 0 1 ; 0.5 0 0 ; 0.5 0 0 ; 0.5 0 0]; %ok

ssbcc(:,:,22) = [-1  1  2 ;  1 -1  1 ; 0 1 0 ; 0 0 0.5 ; 1 0 1 ; 1 1 0.5 ; 1 1 0.5 ; 1 1 0.5];
ssbcc(:,:,23) = [ 1  2  1 ;  1 -1  1 ; 0 1 0 ; 0 0.5 1 ; 1 0 1 ; 1 0.5 0 ; 1 0.5 0 ; 1 0.5 0]; %ok
ssbcc(:,:,24) = [-2 -1  1 ;  1 -1  1 ; 0 1 0 ; 0.5 0 0 ; 1 0 1 ; 0.5 1 1 ; 0.5 1 1 ; 0.5 1 1]; %ok
i112 = 13;
f112 = 24;

% Mode 3,   plane direction, Define four points in the plane 
ssbcc(:,:,25) = [-1 -2  3 ;  1  1  1 ; 0 0 0 ; 1 0 1/3 ; 1 1 1 ; 0 1 2/3 ; 0 1 2/3 ; 0 1 2/3]; %ok
ssbcc(:,:,26) = [-2 -1  3 ;  1  1  1 ; 0 0 0 ; 1 0 2/3 ; 1 1 1 ; 0 1 1/3 ; 0 1 1/3 ; 0 1 1/3]; %ok
ssbcc(:,:,27) = [ 2 -3  1 ;  1  1  1 ; 0 0 0 ; 1 2/3 0 ; 1 1 1 ; 0 1/3 1 ; 0 1/3 1 ; 0 1/3 1]; %ok
ssbcc(:,:,28) = [ 1 -3  2 ;  1  1  1 ; 0 0 0 ; 1 1/3 0 ; 1 1 1 ; 0 2/3 1 ; 0 2/3 1 ; 0 2/3 1]; %ok
ssbcc(:,:,29) = [-3  1  2 ;  1  1  1 ; 0 0 0 ; 2/3 0 1 ; 1 1 1 ; 1/3 1 0 ; 1/3 1 0 ; 1/3 1 0]; %ok
ssbcc(:,:,30) = [-3  2  1 ;  1  1  1 ; 0 0 0 ; 1/3 0 1 ; 1 1 1 ; 2/3 1 0 ; 2/3 1 0 ; 2/3 1 0]; %ok

ssbcc(:,:,31) = [ 1 -2  3 ; -1  1  1 ; 1 0 0 ; 1 1 2/3 ; 0 1 1 ; 0 0 1/3 ; 0 0 1/3 ; 0 0 1/3]; %ok
ssbcc(:,:,32) = [ 2 -1  3 ; -1  1  1 ; 1 0 0 ; 1 1 1/3 ; 0 1 1 ; 0 0 2/3 ; 0 0 2/3 ; 0 0 2/3]; %ok
ssbcc(:,:,33) = [-2 -3  1 ; -1  1  1 ; 1 0 0 ; 1 1/3 1 ; 0 1 1 ; 0 2/3 0 ; 0 2/3 0 ; 0 2/3 0]; %ok
ssbcc(:,:,34) = [-1 -3  2 ; -1  1  1 ; 1 0 0 ; 1 2/3 1 ; 0 1 1 ; 0 1/3 0 ; 0 1/3 0 ; 0 1/3 0]; %ok
ssbcc(:,:,35) = [ 3  1  2 ; -1  1  1 ; 1 0 0 ; 2/3 1 0 ; 0 1 1 ; 1/3 0 1 ; 1/3 0 1 ; 1/3 0 1]; %ok  
ssbcc(:,:,36) = [ 3  2  1 ; -1  1  1 ; 1 0 0 ; 1/3 1 0 ; 0 1 1 ; 2/3 0 1 ; 2/3 0 1 ; 2/3 0 1]; %ok  

ssbcc(:,:,37) = [ 1  2  3 ; -1 -1  1 ; 1 1 0 ; 0 1 1/3 ; 0 0 1 ; 1 0 2/3 ; 1 0 2/3 ; 1 0 2/3]; %ok
ssbcc(:,:,38) = [ 2  1  3 ; -1 -1  1 ; 1 1 0 ; 0 1 2/3 ; 0 0 1 ; 1 0 1/3 ; 1 0 1/3 ; 1 0 1/3]; %ok
ssbcc(:,:,39) = [-2  3  1 ; -1 -1  1 ; 1 1 0 ; 0 1/3 0 ; 0 0 1 ; 1 2/3 1 ; 1 2/3 1 ; 1 2/3 1]; %ok
ssbcc(:,:,40) = [-1  3  2 ; -1 -1  1 ; 1 1 0 ; 0 2/3 0 ; 0 0 1 ; 1 1/3 1 ; 1 1/3 1 ; 1 1/3 1]; %ok
ssbcc(:,:,41) = [ 3 -1  2 ; -1 -1  1 ; 1 1 0 ; 1/3 1 1 ; 0 0 1 ; 2/3 0 0 ; 2/3 0 0 ; 2/3 0 0]; %ok
ssbcc(:,:,42) = [ 3 -2  1 ; -1 -1  1 ; 1 1 0 ; 2/3 1 1 ; 0 0 1 ; 1/3 0 0 ; 1/3 0 0 ; 1/3 0 0]; %ok

ssbcc(:,:,43) = [-1  2  3 ;  1 -1  1 ; 0 1 0 ; 0 0 2/3 ; 1 0 1 ; 1 1 1/3 ; 1 1 1/3 ; 1 1 1/3]; %ok
ssbcc(:,:,44) = [-2  1  3 ;  1 -1  1 ; 0 1 0 ; 0 0 1/3 ; 1 0 1 ; 1 1 2/3 ; 1 1 2/3 ; 1 1 2/3]; %ok
ssbcc(:,:,45) = [ 2  3  1 ;  1 -1  1 ; 0 1 0 ; 0 2/3 1 ; 1 0 1 ; 1 1/3 0 ; 1 1/3 0 ; 1 1/3 0]; %ok
ssbcc(:,:,46) = [ 1  3  2 ;  1 -1  1 ; 0 1 0 ; 0 1/3 1 ; 1 0 1 ; 1 2/3 0 ; 1 2/3 0 ; 1 2/3 0]; %ok
ssbcc(:,:,47) = [-3 -1  2 ;  1 -1  1 ; 0 1 0 ; 1/3 0 0 ; 1 0 1 ; 2/3 1 1 ; 2/3 1 1 ; 2/3 1 1]; %ok
ssbcc(:,:,48) = [-3 -2  1 ;  1 -1  1 ; 0 1 0 ; 2/3 0 0 ; 1 0 1 ; 1/3 1 1 ; 1/3 1 1 ; 1/3 1 1]; %ok
i123 = 25;
f123 = 48;

% FCC
% Mode 1,  plane direction, Define four points in the plane       For TiAl, 1,4,7,10 are ordinary dislocations
ssfcc(:,:,1) =  [ 1  1  1 ; -1  1  0 ; 1 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ; 0 0 1 ; 0 0 1];  %ok  others are superdislocations
ssfcc(:,:,2) =  [ 1  1  1 ;  1  0 -1 ; 0 0 1 ; 0 0 1 ; 1 0 0 ; 0 1 0 ; 0 1 0 ; 0 1 0];  %ok   13, 16, 19, 22 are true twins
ssfcc(:,:,3) =  [ 1  1  1 ;  0 -1  1 ; 0 1 0 ; 0 1 0 ; 0 0 1 ; 1 0 0 ; 1 0 0 ; 1 0 0];  %ok    others are pseudo twins/superdis

ssfcc(:,:,4) =  [-1  1  1 ; -1 -1  0 ; 1 1 0 ; 1 1 0 ; 0 0 0 ; 1 0 1 ; 1 0 1 ; 1 0 1];  %ok
ssfcc(:,:,5) =  [-1  1  1 ;  1  0  1 ; 0 0 0 ; 0 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 0 ; 1 1 0];  %ok
ssfcc(:,:,6) =  [-1  1  1 ;  0  1 -1 ; 1 0 1 ; 1 0 1 ; 1 1 0 ; 0 0 0 ; 0 0 0 ; 0 0 0];  %ok

ssfcc(:,:,7) =  [-1 -1  1 ;  1 -1  0 ; 0 1 0 ; 0 1 0 ; 1 0 0 ; 1 1 1 ; 1 1 1 ; 1 1 1];  %ok
ssfcc(:,:,8) =  [-1 -1  1 ; -1  0 -1 ; 1 1 1 ; 1 1 1 ; 0 1 0 ; 1 0 0 ; 1 0 0 ; 1 0 0];  %ok
ssfcc(:,:,9) =  [-1 -1  1 ;  0  1  1 ; 1 0 0 ; 1 0 0 ; 1 1 1 ; 0 1 0 ; 0 1 0 ; 0 1 0];  %ok

ssfcc(:,:,10) = [ 1 -1  1 ;  1  1  0 ; 0 0 0 ; 0 0 0 ; 1 1 0 ; 0 1 1 ; 0 1 1 ; 0 1 1];  %ok
ssfcc(:,:,11) = [ 1 -1  1 ; -1  0  1 ; 1 1 0 ; 1 1 0 ; 0 1 1 ; 0 0 0 ; 0 0 0 ; 0 0 0];  %ok
ssfcc(:,:,12) = [ 1 -1  1 ;  0 -1 -1 ; 0 1 1 ; 0 1 1 ; 0 0 0 ; 1 1 0 ; 1 1 0 ; 1 1 0];  %ok

% Mode 2 twins    plane    direction, Define four points in the plane 
ssfcc(:,:,13) = [ 1  1  1 ;  1  1 -2 ; 0 0 1 ; 1 0 0 ;  0.5  0.5  0 ; 0 1 0 ; 0 1 0 ; 0 1 0];  %ok
ssfcc(:,:,14) = [ 1  1  1 ;  1 -2  1 ; 0 1 0 ; 0 0 1 ;  0.5  0  0.5 ; 1 0 0 ; 1 0 0 ; 1 0 0];  %ok
ssfcc(:,:,15) = [ 1  1  1 ; -2  1  1 ; 1 0 0 ; 0 1 0 ;  0  0.5  0.5 ; 0 0 1 ; 0 0 1 ; 0 0 1];  %ok

ssfcc(:,:,16) = [-1  1  1 ; -1  1 -2 ; 1 0 1 ; 1 1 0 ;  0.5  0.5  0 ; 0 0 0 ; 0 0 0 ; 0 0 0];  %ok
ssfcc(:,:,17) = [-1  1  1 ; -1 -2  1 ; 1 1 0 ; 0 0 0 ;  0.5  0  0.5 ; 1 0 1 ; 1 0 1 ; 1 0 1];  %ok
ssfcc(:,:,18) = [-1  1  1 ;  2  1  1 ; 0 0 0 ; 1 0 1 ;  1  0.5  0.5 ; 1 1 0 ; 1 1 0 ; 1 1 0];  %ok

ssfcc(:,:,19) = [-1 -1  1 ; -1 -1 -2 ; 1 1 1 ; 0 1 0 ;  0.5  0.5  0 ; 1 0 0 ; 1 0 0 ; 1 0 0];  %ok
ssfcc(:,:,20) = [-1 -1  1 ; -1  2  1 ; 1 0 0 ; 1 1 1 ;  0.5  1  0.5 ; 0 1 0 ; 0 1 0 ; 0 1 0];  %ok
ssfcc(:,:,21) = [-1 -1  1 ;  2 -1  1 ; 0 1 0 ; 1 0 0 ;  1  0.5  0.5 ; 1 1 1 ; 1 1 1 ; 1 1 1];  %ok

ssfcc(:,:,22) = [ 1 -1  1 ;  1 -1 -2 ; 0 1 1 ; 0 0 0 ;  0.5  0.5  0 ; 1 1 0 ; 1 1 0 ; 1 1 0];  %ok
ssfcc(:,:,23) = [ 1 -1  1 ;  1  2  1 ; 0 0 0 ; 1 1 0 ;  0.5  1  0.5 ; 0 1 1 ; 0 1 1 ; 0 1 1];  %ok
ssfcc(:,:,24) = [ 1 -1  1 ;  2  1 -1 ; 1 1 0 ; 0 1 1 ;  0  0.5  0.5 ; 0 0 0 ; 0 0 0 ; 0 0 0];  %ok

% BCT (Sn)
% Mode 1,      plane direction, Define four points in the plane 
ssbct(:,:,1) = [1,0,0 ; 0,0,1; 0 1 0 ; 0 1 1 ; 0 1 1 ; 0 0 1 ; 0 0 1; 0 1 0];
ssbct(:,:,2) = [0,1,0 ; 0,0,1; 0 0 0 ; 0 0 1 ; 0 0 1 ; 1 0 1 ; 1 0 0; 0 0 0];
% Mode 2,      first point of 4 is base for Burger's vector and plane normal
ssbct(:,:,3) = [1,1,0 ; 0,0,1 ; 1 0 0 ; 1 0 1 ; 1 0 1 ; 0 1 1 ; 0 1 0 ; 1 0 0];
ssbct(:,:,4) = [1,-1,0 ; 0,0,1 ; 0 0 0 ; 1 1 0 ; 1 1 1 ; 0 0 1 ; 0 0 0 ; 0 0 0];
% Mode 3
ssbct(:,:,5) = [1,0,0 ; 0,1,0 ; 0 0 0 ; 0 0 0 ; 0 1 0 ; 0 1 1 ; 0 0 1; 0 0 0];
ssbct(:,:,6) = [0,1,0 ; 1,0,0 ; 0 0 0 ; 0 0 0 ; 1 0 0 ; 1 0 1 ; 0 0 1; 0 0 0];
% Mode 4
ssbct(:,:,7) = [1,1,0 ; 1,-1,1 ; 0 1 0 ; 1 0 0 ; 1 0 1 ; 0 1 1; 0 1 0 ; 0 1 0 ];
ssbct(:,:,8) = [1,1,0 ; -1,1,1 ; 1 0 0 ; 0 1 0 ; 0 1 1 ; 1 0 1 ; 1 0 0 ; 1 0 0];
ssbct(:,:,9) = [1,-1,0 ;  1,1,1 ;  0 0 0 ; 1 1 0 ; 1 1 1 ; 0 0 1; 0 0 0; 0 0 0];
ssbct(:,:,10) = [1,-1,0 ; -1,-1,1 ; 1 1 0 ; 1 1 1 ; 0 0 1; 0 0 0 ; 1 1 0 ; 1 1 0];
% Mode 5
ssbct(:,:,11) = [1,1,0 ; 1,-1,0 ; 0 1 0 ; 1 0 0 ; 1 0 1 ; 0 1 1 ; 0 1 0 ; 0 1 0];
ssbct(:,:,12) = [1,-1,0 ; 1,1,0 ; 0 0 0 ; 1 1 0 ; 1 1 1 ; 0 0 1 ; 0 0 0 ; 0 0 0];
% Mode 6
ssbct(:,:,13) = [0,1,0 ; 1,0,1  ; 0 0 0 ; 0 0 1 ; 1 0 1 ; 1 0 0; 0 0 0; 0 0 0];
ssbct(:,:,14) = [0,1,0 ; 1,0,-1 ; 0 0 1 ; 1 0 1 ; 1 0 0 ; 0 0 0; 0 0 1; 0 0 1];
ssbct(:,:,15) = [1,0,0 ; 0,1,1  ; 0 0 0 ; 0 1 0 ; 0 1 1 ; 0 0 1; 0 0 0; 0 0 0];
ssbct(:,:,16) = [1,0,0 ; 0,1,-1 ; 0 0 1 ; 0 0 0 ; 0 1 0 ; 0 1 1; 0 0 1; 0 0 1];
% Mode 7
ssbct(:,:,17) = [0,0,1 ; 1,0,0 ; 0 0 0 ; 1 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0; 0 0 0];
ssbct(:,:,18) = [0,0,1 ; 0,1,0 ; 0 0 0 ; 0 1 0 ; 0 1 0 ; 1 1 0 ; 1 0 0; 0 0 0];
% Mode 8
ssbct(:,:,19) = [0,0,1 ; 1,1,0  ; 0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 0 0 ; 0 0 0];
ssbct(:,:,20) = [0,0,1 ; 1,-1,0 ; 0 1 0 ; 0 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 0 1 0];
% Mode 9
ssbct(:,:,21) = [1,0,1 ; 1,0,-1 ; 0 0 1 ; 1 0 0 ; 1 1 0 ; 0 1 1 ; 0 0 1 ; 0 0 1];
ssbct(:,:,22) = [1,0,-1 ; 1,0,1 ; 0 0 0 ; 0 1 0 ; 1 1 1 ; 1 0 1 ; 0 0 0 ; 0 0 0];
ssbct(:,:,23) = [0,1,1 ; 0,1,-1 ; 0 0 1 ; 1 0 1 ; 1 1 0 ; 0 1 0 ; 0 0 1 ; 0 0 1];
ssbct(:,:,24) = [0,1,-1 ; 0,1,1 ; 0 0 0 ; 1 0 0 ; 1 1 1 ; 0 1 1 ; 0 0 0 ; 0 0 0];
% Mode 10,  three points are defined in the plane (one repeated)
ssbct(:,:,25) = [1,2,1 ; -1,0,1 ;  1 0 0 ; 0 .5 0 ; 0 0 1 ; 1 0 0 ; 1 0 0 ; 1 0 0];
ssbct(:,:,26) = [1,-2,1 ; -1,0,1 ; 1 1 0 ; 1 1 0 ; 0 1 1 ; 0 .5 0 ; 1 1 0 ; 1 1 0];
ssbct(:,:,27) = [-1,2,1 ; 1,0,1 ;  0 0 0 ; 0 0 0 ; 1 0 1 ; 1 .5 0 ; 0 0 0 ; 0 0 0];
ssbct(:,:,28) = [-1,-2,1 ; 1,0,1 ; 0 1 0 ; 1 .5 0 ; 1 1 1 ; 0 1 0 ; 0 1 0 ; 0 1 0];
ssbct(:,:,29) = [2,1,1 ; 0,-1,1 ;  0 1 0 ; 0 1 0 ; 0 0 1 ; .5 0 0 ; 0 1 0 ; 0 1 0];
ssbct(:,:,30) = [-2,1,1 ; 0,-1,1 ; 1 1 0 ; .5 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 0 ; 1 1 0];
ssbct(:,:,31) = [-2,-1,1 ; 0,1,1 ; 1 0 0 ; 1 0 0 ; 1 1 1 ; .5 1 0 ; 1 0 0 ; 1 0 0];
ssbct(:,:,32) = [2,-1,1 ; 0,1,1 ;  0 0 0 ; .5 1 0 ; 0 1 1 ; 0 0 0 ; 0 0 0 ; 0 0 0];

% Note that plotting may not work yet for Sn...

mnslp = max(nslp);
ss = zeros(8,3,mnslp,4);

for i=1:1:mnslp   % Change n & m to unit vector,
    if i <= nslphex
        n=[sshex(1,1,i)  (sshex(1,2,i)*2+sshex(1,1,i))/3^.5  sshex(1,4,i)/c_a_hex]; % Plane normal /c_a_hex
        m=[sshex(2,1,i)*1.5  3^.5/2*(sshex(2,2,i)*2+sshex(2,1,i))  sshex(2,4,i)*c_a_hex]; % Slip direction *c_a_hex 
        ss(1,:,i,1) = n/norm(n);  % alpha plane
        ss(2,:,i,1) = m/norm(m);  % alpha direction                               PHASE 1 is HEX
        ss(3,:,i,1) = [3*sshex(3,1,i)/2, (sshex(3,1,i)+2*sshex(3,2,i))*sqrt(3)/2, sshex(3,4,i)*c_a_hex]; % hpoint 1
        ss(4,:,i,1) = [3*sshex(4,1,i)/2, (sshex(4,1,i)+2*sshex(4,2,i))*sqrt(3)/2, sshex(4,4,i)*c_a_hex]; % hpoint 2
        ss(5,:,i,1) = [3*sshex(5,1,i)/2, (sshex(5,1,i)+2*sshex(5,2,i))*sqrt(3)/2, sshex(5,4,i)*c_a_hex]; % hpoint 3
        ss(6,:,i,1) = [3*sshex(6,1,i)/2, (sshex(6,1,i)+2*sshex(6,2,i))*sqrt(3)/2, sshex(6,4,i)*c_a_hex]; % hpoint 4
        ss(7,:,i,1) = [3*sshex(7,1,i)/2, (sshex(7,1,i)+2*sshex(7,2,i))*sqrt(3)/2, sshex(7,4,i)*c_a_hex]; % hpoint 5
        ss(8,:,i,1) = [3*sshex(8,1,i)/2, (sshex(8,1,i)+2*sshex(8,2,i))*sqrt(3)/2, sshex(8,4,i)*c_a_hex]; % hpoint 6
    end

    if i <= nslpbcc
        n = [ssbcc(1,1,i),ssbcc(1,2,i),ssbcc(1,3,i)/c_a_bcc];   % slightly tetragonal has c/a <> 1.0
        m = [ssbcc(2,1,i),ssbcc(2,2,i),ssbcc(2,3,i)*c_a_bcc]; 
        ss(1,:,i,2) = n/norm(n);  % bcc plane                                     PHASE 2 is BCC
        ss(2,:,i,2) = m/norm(m);  % bcc direction
        ss(3,:,i,2) = [ssbcc(3,1,i),ssbcc(3,2,i),ssbcc(3,3,i)*c_a_bcc]; % point 1
        ss(4,:,i,2) = [ssbcc(4,1,i),ssbcc(4,2,i),ssbcc(4,3,i)*c_a_bcc]; % point 2
        ss(5,:,i,2) = [ssbcc(5,1,i),ssbcc(5,2,i),ssbcc(5,3,i)*c_a_bcc]; % point 3
        ss(6,:,i,2) = [ssbcc(6,1,i),ssbcc(6,2,i),ssbcc(6,3,i)*c_a_bcc]; % point 4
        ss(7,:,i,2) = [ssbcc(7,1,i),ssbcc(7,2,i),ssbcc(7,3,i)*c_a_bcc]; % point 5
        ss(8,:,i,2) = [ssbcc(8,1,i),ssbcc(8,2,i),ssbcc(8,3,i)*c_a_bcc]; % point 6
    end
    
    if i <= nslpfcc
        n = [ssfcc(1,1,i),ssfcc(1,2,i),ssfcc(1,3,i)/c_a_fcc];   % slightly tetragonal has c/a <> 1.0
        m = [ssfcc(2,1,i),ssfcc(2,2,i),ssfcc(2,3,i)*c_a_fcc]; 
        ss(1,:,i,3) = n/norm(n);  % fcc plane                                     PHASE 3 is FCC
        ss(2,:,i,3) = m/norm(m);  % fcc direction
        ss(3,:,i,3) = [ssfcc(3,1,i),ssfcc(3,2,i),ssfcc(3,3,i)*c_a_fcc]; % point 1
        ss(4,:,i,3) = [ssfcc(4,1,i),ssfcc(4,2,i),ssfcc(4,3,i)*c_a_fcc]; % point 2
        ss(5,:,i,3) = [ssfcc(5,1,i),ssfcc(5,2,i),ssfcc(5,3,i)*c_a_fcc]; % point 3
        ss(6,:,i,3) = [ssfcc(6,1,i),ssfcc(6,2,i),ssfcc(6,3,i)*c_a_fcc]; % point 4
        ss(7,:,i,3) = [ssfcc(7,1,i),ssfcc(7,2,i),ssfcc(7,3,i)*c_a_fcc]; % point 5
        ss(8,:,i,3) = [ssfcc(8,1,i),ssfcc(8,2,i),ssfcc(8,3,i)*c_a_fcc]; % point 6
    end
    
    if i <= nslpbct
        n = [ssbct(1,1,i),ssbct(1,2,i),ssbct(1,3,i)/c_a_bct]; 
        m = [ssbct(2,1,i),ssbct(2,2,i),ssbct(2,3,i)*c_a_bct]; 
        ss(1,:,i,4) = n/norm(n);  % bct plane                                     PHASE 4 is BCT
        ss(2,:,i,4) = m/norm(m);  % bct direction
        ss(3,:,i,4) = [ssbct(3,1,i),ssbct(3,2,i),ssbct(3,3,i)*c_a_bct]; % point 1
        ss(4,:,i,4) = [ssbct(4,1,i),ssbct(4,2,i),ssbct(4,3,i)*c_a_bct]; % point 2
        ss(5,:,i,4) = [ssbct(5,1,i),ssbct(5,2,i),ssbct(5,3,i)*c_a_bct]; % point 3
        ss(6,:,i,4) = [ssbct(6,1,i),ssbct(6,2,i),ssbct(6,3,i)*c_a_bct]; % point 4
        ss(7,:,i,4) = [ssbct(7,1,i),ssbct(7,2,i),ssbct(7,3,i)*c_a_bct]; % point 5
        ss(8,:,i,4) = [ssbct(8,1,i),ssbct(8,2,i),ssbct(8,3,i)*c_a_bct]; % point 6
    end

%   translator from old cell version to new ss() structure    
%     n =[slpsys{2,iss}(1,1),slpsys{2,iss}(1,2),slpsys{2,iss}(1,3)/c_a_bct]; % Plane normal  /c_a_bct
%     m =[slpsys{2,iss}(2,1),slpsys{2,iss}(2,2),slpsys{2,iss}(2,3)*c_a_bct]; % Slip direction  *c_a_bct
%     unit_n = n/norm(n);
%     unit_m = m/norm(m);
%     p1=[slpsys{2,iss}(3,1),slpsys{2,iss}(3,2),slpsys{2,iss}(3,3)*c_a_bct]; % point 1
%     p2=[slpsys{2,iss}(4,1),slpsys{2,iss}(4,2),slpsys{2,iss}(4,3)*c_a_bct]; % point 2
%     p3=[slpsys{2,iss}(5,1),slpsys{2,iss}(5,2),slpsys{2,iss}(5,3)*c_a_bct]; % point 3
%     p4=[slpsys{2,iss}(6,1),slpsys{2,iss}(6,2),slpsys{2,iss}(6,3)*c_a_bct]; % point 4
end


% %  Loop for grains to establish slip conditions for each grain


Fgrcen = zeros(int16(dIDgr(1,1)*1.1),18);  % this sets up an array for grain information 
EY = zeros(int16(dIDgr(1,1)*1.1),3);
Sfplbv = zeros(mnslp+1,30);
sortmv = zeros(mnslp+1,30,int16(dIDgr(1,1)*1.1));
grcen(:,1) = -1;    % that is a little bigger that needed because some grain numbers are skipped, 
                    % and are thus marked with -1.  Grains are processed by grain number, not array location
grmax = 0;    ngcount = 0;   fprintf('Numbers and vectors computed for Grain # ');
for ng=1:1:dIDgr(1,1); %    generalized Schmid factor calculation loop for each grain ng
    if ng>ngcount+dIDgr/10;
        ngcount=ngcount+dIDgr/10;
        fprintf(' %d ',ng);
    end
    ig = IDgr(ng,1);
    if ig > grmax
        grmax = ig;
    end           %  phase(1)  grain center(2,3)   eulers(4:6)
    if ig > 0
        grcen(ig,1:6) = [IDgr(ng,10) IDgr(ng,5:6) IDgr(ng,2:4)];
        phid = grcen(ig,4:6);   % phid is Euler phi angles in degrees

        ph = grcen(ig,1);   % phase ID set    

        if hkl == 1
            phid(1) = phid(1) + 0 ; % or +90 to convert hkl to TSL software default   
            if phid(1)>360          % or +180 to modify TSL Euler angle coordinate system to have X down and Y right;
                phid(1) = phid(1) - 360;
            elseif phid(1) < 0
                phid(1) = phid(1) + 360;
            end
        end

        g1=[cosd(phid(1)),sind(phid(1)),0; -sind(phid(1)),cosd(phid(1)),0; 0,0,1];
        g2=[1,0,0; 0,cosd(phid(2)),sind(phid(2)); 0,-sind(phid(2)),cosd(phid(2))];
        g3=[cosd(phid(3)),sind(phid(3)),0; -sind(phid(3)),cosd(phid(3)),0; 0,0,1];
        g=g3*g2*g1;
        if nsten == 1
            sigma_n(:,:,ig) = sigma_n(:,:,1);
            sigma_v(:,ig) = sigma_v(:,1);
        end
        gsgT = g*sigma_n(:,:,ig)*g'; %rotated stress tensor
        grcen(ig,7:9) = [0 0 1]*g; %c-axis direction
        grcen(ig,10:18) = [g(1,:) g(2,:) g(3,:)];  % Orientation matrix is stored

    % calculate elastic modulus to find compliance mismatch in three principal directions  (from Nye textbook on Anisotropy)   
        if ph == 5  % needs a different structure for the sij matrix - needs more terms, not working in this version. 
            % E(orth) = ( s11*x^4 + 2s12*x^2*y^2 + 2s13*x^2*y^2 + s22*y^4 + 2s23*y^2*z^2 + s33*z^4 + s44*y^2*z^2 + s55*x^2*z^2 + s66*y^2*z^2 )^-1*100
        elseif ph == 4
            % E(tetr) = ( s11 * (x^4 + y^4) + s33 * z^4 + (s66 + 2*s12)*x^2*y^2 + (2*s13 + s44) * z^2 * (y^2* + x^2) )^-1*100
            for i = 1:1:3   %  This gives the modulus in the x (1), y (2), z (3) directions (as looped by i)
    %         CTExyz(ng,i) = (CTEm(1,i)^2+CTEm(2,i)^2+CTEm(3,i)^2)^.5;   %   coefficient of thermal expansion anisotropy
            e1 = sij(ph,1)*(g(1,i)^4 + g(2,i)^4) + sij(ph,4)*g(3,i)^4;
            e2 = (2*sij(ph,2) + sij(ph,6))*(g(1,i)^2 * g(2,i)^2);
            e3 = (2*sij(ph,3) + sij(ph,5))*g(3,i)^2 * (g(1,i)^2 + g(2,i)^2);
            EY(ig,i)=1./(e1+e2+e3);   % NOTE: slip system and plane information not installed in this version for Sn or TiAl
            end
        elseif ph == 2 || ph == 3
            % E(cube) = ( s11 - 2*(s11 - s12 - s44/2)*(x^2*y^2 + y^2*z^2 + z^2*x^2) )^-1
            for i = 1:1:3
            EY(ig,i) = 1./(sij(ph,1) - 2*(sij(ph,1) - sij(ph,2) - sij(ph,5)/2) *...
                (g(1,i)^2*g(2,i)^2 + g(2,i)^2*g(3,i)^2 + g(3,i)^2*g(1,i)^2) );
            end
        elseif ph == 1
            % E(hex)  = ( s11 * (1-z^2)^2 + s33 * z^4 + (s44 + 2*s13)*(1-z^2)*z^2 )^-1
            for i = 1:1:3
            EY(ig,i) = 1./(sij(ph,1) * (1-g(3,i)^2)^2 + sij(ph,4) * g(3,i)^4 +...
                (sij(ph,5) + 2*sij(ph,3))*(1-g(3,i)^2)*g(3,i)^2 );
            end
        end

        for j=1:1:nslp(ph)  %  direction          plane    Sfplbv means Schmid factor, plane and Burgers vector (and points on plane)
            Sfplbv(j,1) = j;   %  n * sigma * m
            Sfplbv(j,2) = ss(2,:,j,ph)*gsgT*ss(1,:,j,ph)';  %  generalized Schmid factor
            if ph == 1 && j>24 && Sfplbv(j,2)<0 
                Sfplbv(j,2) = 0.001*Sfplbv(j,2) ;    % this is to prevent anti-twin shears from being seriously considered later
            end
            Sfplbv(j,3) = abs(Sfplbv(j,2));   % abs(generalized schmid factor)
            Sfplbv(j,4:6) = g'*ss(1,:,j,ph)'; % plane normal in lab coords
            Sfplbv(j,7:9) = g'*ss(2,:,j,ph)'; % bv direction in lab coords
            Sfplbv(j,10:12) = cross(Sfplbv(j,4:6),[0,0,1]); % plane trace
            for k = 1:1:6
                is = 3*k+10;
                ie = is+2;
                Sfplbv(j,is:ie) = g'*ss(k+2,:,j,ph)'; % plane plotting vectors from origin to points in cell, in lab coords
            end
        end                        %useful plotting for hexahedral tetragonal unit cell vectors that sort to bottom row
        Sfplbv(mnslp+1,:) = [ph 1 -1 [1 0 0]*g [0 1 0]*g [0 0 1*c_a(ph)]*g 0 0 0  0 0 0  0 0 0  0 0 0  0 0 0  0 0 0]; % don't change g to g' here! otherwise it may make incorrect cubic prisms
        sortmv(:,:,ig)  = sortrows(Sfplbv,-3);  % Sort slip systems by Schimd factor
    end  %  ig > 0 check
end  % ng loop
fprintf(' %d\n ', ng);


% %  Now start processing by grain boundary...


gbnorm = zeros(dRBdy(1,1),3); 
gbtrac = zeros(dRBdy(1,1),3);
damp   = zeros(22,dRBdy(1,1));
mpr    = zeros(mnslp+1,mnslp+1,dRBdy(1,1));
nbcount = 0;  fprintf('Computing grain boundary parameters for GB # ');
jk4max = 1; top3d3 = 0; top3d2 = 0; top3d1 = 0; top3d0 = 0; 
for gbnum = 1:1:dRBdy(1,1);  %gbnum is grain boundary number, will calculate m' and other damage parameters
    if gbnum>nbcount+dRBdy/10;
        nbcount=nbcount+dRBdy/10;
        fprintf(' %d ',gbnum);     %, jk4max
    end
    
    grA = RBdy(gbnum,13);    grB = RBdy(gbnum,14);
    jk = 1;   % counter for the number of m' calculations made where Schmid factors are > low tolerance
    jk4 = 1;  % counter for the number of m' calculations made where Schmid factors are > high tolerance 
    mpmax = 0;    mploc = 0;
    dpsum = 0;    dpsum4 = 0;
    mpsum = 0;    mpsum4 = 0;
    damage = 0;   damage4 = 0;
    
    mpr(1,1,gbnum) = grA+grB/1000; % m-prime table label for grain numbers in mpr()

    if strcmp(num2str(sigma_v(:,ig)),num2str([1 0 0]')) % stress axis || [100] (X)      
        EgrA = EY(grA,1); EgrB = EY(grB,1);        
    elseif strcmp(num2str(sigma_v(:,ig)),num2str([0 1 0]')) % stress axis || [010] (Y)
        EgrA = EY(grA,2); EgrB = EY(grB,2);        
    elseif strcmp(num2str(sigma_v(:,ig)),num2str([0 0 1]'))  % stress axis || [001] (Z)  
        EgrA = EY(grA,3); EgrB = EY(grB,3);                
    elseif trace(sigma_n(:,:,ig)) == 0    % Crude estimate of shear effects follows, may not be meaningful
            EgrA = EY(grA,3)*abs(sigma_n(1,2,ig)) + EY(grA,2)*abs(sigma_n(1,3,ig)) + EY(grA,1)*abs(sigma_n(2,3,ig)); 
            EgrB = EY(grB,3)*abs(sigma_n(1,2,ig)) + EY(grB,2)*abs(sigma_n(1,3,ig)) + EY(grB,1)*abs(sigma_n(2,3,ig));
    else
        fprintf('Can''t calculate modulus for this stress state\r');
        EgrA = 1; EgrB = 1  %pause
    end
        Eratio = min(EgrA,EgrB)/max(EgrA,EgrB);        % always use Emin/Emax!
    
    F1A = zeros(1,mnslp); % F1 FIP, Simkin et al. 2003 for grain A 
%     F14A = 0; % F1 FIP w/ restriction on Schmid factor value for grain A 
    F1B = zeros(1,mnslp); % F1 FIP for grain B
%     F14B = 0; % F1 FIP w/ restriction on Schmid factor value for grain B
    F1 = 0; % F1 for grainA/grainB
    F14A = zeros(1,mnslp); % F14 for grainA (with restriction on Schmid factor value)
    F14B = zeros(1,mnslp);
    
    normdp4 = .75;  % These values indicate instances where values 
    normp4 = .55;   % of dm' or m' are too low to take seriously
%    RBdy(gbnum,1:6) = (180/pi).*RBdy(gbnum,1:6);
    gbnorm(gbnum,:) = [cosd(RBdy(gbnum,8)) sind(RBdy(gbnum,8)) 0];
    gbtrac(gbnum,:) = sigma_n(:,:,grA)*gbnorm(gbnum,:)';
    
    for k = 1:1:mnslp                        % Build table of m' values for each grain pair
        schmA = sortmv(k,2,grA);  % Use the correctly signed version of Schmid factor in position 2
        mpr(k+1,1,gbnum) = sortmv(k,1,grA) + abs(schmA) ;  % schmid factor header  grA is down (k goes with A)
        for j = 1:1:mnslp  % j goes across (with grain B)
            schmB = sortmv(j,2,grB);  % Use the correctly signed version of Schmid factor in position 2
            mpr(1,j+1,gbnum) = sortmv(j,1,grB) + abs(schmB) ;  % grA goes down, grB across, (j goes with B)
                if abs(schmA) > schmid_tolL && abs(schmB) > schmid_tolL
                    mp = sortmv(j,4:6,grB)*sortmv(k,4:6,grA)'; %  plane        This is a dot product
                    mb = sortmv(j,7:9,grB)*sortmv(k,7:9,grA)'*sign(schmA*schmB); %  direction    This is a dot product
                    mprime = mp*mb;                        %  if the b directions are pointing similarly mb must be positive
                    mpr(k+1,j+1,gbnum) = mprime;   % m' table is filled in, next do some things with this value...
                                                % a damage parameter is found assuming m' = 0.8 is worst case                                                                                       
                                           % schmA                 
                    F1A(k) = F1A(k)+abs(sortmv(k,3,grA)*dot(sortmv(k,7:9,grA)',sigma_v(:,ig))*dot(sortmv(k,7:9,grA)',sortmv(j,7:9,grB)')); % F1 for grain A
                    F1B(j) = F1B(j)+abs(sortmv(j,3,grB)*dot(sortmv(j,7:9,grB)',sigma_v(:,ig))*dot(sortmv(j,7:9,grB)',sortmv(k,7:9,grA)')); % F1 for grain B

                    if mprime > 0.6  % assumes no slip transfer (or damage) occurs if m' < 0.6 or when m' = 1
                        dampar = mprime + 0.1;   % pushes m' up by 0.1, so that 1 is worst case.
                        if dampar > 1                   
                            dampar = 2-dampar;          % assume dampar < 1 means less likely to generate damage 
                        end
                        dpsum = dpsum + dampar;         % damage parameter only, for all slip systems where m' > 0.8
                        mpsum = mpsum + abs(mprime);         % m'  only, for all slip systems where m' > 0.8
                        damage = damage + dampar*max(abs(schmA),abs(schmB));   % damage parameter modified by schmid factor

                        if (abs(schmA) > schmid_tolH && abs(schmB) > schmid_tolH) || abs(schmA) < 0.001 || abs(schmB) < 0.001
                            %  assumes that slip transfer happens only if m' > 0.35 and both have high schmid factor and capturing effect of anti-twin
                            himp4(jk4,:,gbnum) = [k j abs(mprime), mprime]; % location and values of high m' values for grain pair; k(row) is from grain c13 j(column) is from grain c14 in RBdy
                            dpsum4 = dpsum4 + dampar;   %  running sum of slip system interactions that have m > 0.4
                            mpsum4 = mpsum4 + abs(mprime);   %  running sum of m' for slip system interactions that have m > 0.4
                            damage4 = damage4 + dampar*max(abs(schmA),abs(schmB));  % damage parameter modified by schmid factor
                            if jk4 > jk4max
                                jk4max = jk4;
                            end
                            jk4 = jk4 + 1;
                            F14A(k) = F14A(k)+abs(sortmv(k,3,grA)*dot(sortmv(k,7:9,grA)',sigma_v(:,ig))*dot(sortmv(k,7:9,grA)',sortmv(j,7:9,grB)')); % F14 for grain A
                            F14B(j) = F14B(j)+abs(sortmv(j,3,grB)*dot(sortmv(j,7:9,grB)',sigma_v(:,ig))*dot(sortmv(j,7:9,grB)',sortmv(k,7:9,grA)')); % F14 for grain B                      
                        end
                        if mpmax < abs(mprime)  % record maximum m' value and its location in matrix
                            mpmax = abs(mprime);
                        end
                        jk = jk + 1;
                    end  % of m' > 0.6 if statement
                end  % of if statement for values with schmid factors > schmid_tolL
        end % slip system j loop for each B grain in pair
    end   % slip system k loop for each A grain in pair
%     Variables evaluated above in the loops:     dpsum     mpsum     jk4     dpsum4     mpsum4     F1     F14
    gbodam = damage*abs(norm(gbtrac(gbnum,:))) ;  % damage parameter modified by apparent GB inclination
    gbodam4 = damage4*abs(norm(gbtrac(gbnum,:))) ; % damage parameter for high schmid modified by GB inclination
    caxmis = (1-(grcen(grA,7:9)*grcen(grB,7:9)')^2)^.5;  % misorientation of c-axes (not meaningful for cubic)
    normdp = dpsum/(jk-1);  % average value of dm'
    normp = mpsum/(jk-1);  % average value of m'
    F1Asort = sort(F1A,'descend'); F1Bsort = sort(F1B,'descend');
    F14Asort = sort(F14A,'descend'); F14Bsort = sort(F14B,'descend');    
    maxF1 = max(F1Asort(1),F1Bsort(1));
    maxF14 = max(F14Asort(1),F14Bsort(1));
    if jk4> 1 
        normdp4 = dpsum4/(jk4-1);  %average value of damage for slip systems with m > 0.4
        normp4 = mpsum4/(jk4-1);  %average value of m' for slip systems with m > 0.4
    end
    
    jk6 = 0; jktop36 = 0;  % counters for number of m' values to average later.
    top6mpn = 0; top3mpn = 0; top3mp6 = 0;
    
    Spair = zeros(5:5); flag = 1;% Spairmp = zeros(5:5);
    for k = 1:1:5        % Rather than finding all m' values for m > tolH, look only for top 3 or top 6
        for j = 1:1:5    % find sum of schmid factors for each element of mpr array and put in Spair()
            Spair(k,j) = mpr(k+1,1,gbnum)-fix(mpr(k+1,1,gbnum)) + mpr(1,j+1,gbnum)-fix(mpr(1,j+1,gbnum));
%            Spairmp(k,j) = Spair(k,j);
        end
    end
    while (jk6 < 6 || jktop36 < 3) && flag > 0
        maxij = [0 0 0];
        for k = 1:1:5
            for j = 1:1:5
                if Spair(k,j)>maxij(1)  %finds largest Spair value in 5x5 Spair()
                    maxij = [Spair(k,j) k j]; %identifies location k j of sum of Schmid factors
                end
            end
        end
        flag = maxij(1);
        mpchk = abs(mpr(maxij(2)+1,maxij(3)+1,gbnum));  % puts (next) highest m' into mpchk
        if jk6 < 6 && flag > 0
            top6mpn = top6mpn + mpchk;
            Spair(maxij(2),maxij(3)) = -1;  % Now that highest Spair value is found and recorded, wipe it out
            if jk6 == 3;
                top3mpn = top6mpn;  %  capture top 3 values in this variable
            end
            jk6 = jk6 + 1;
        end
        if jktop36 < 3 && flag > 0
            Spair(maxij(2),maxij(3)) = Spair(maxij(2),maxij(3)) -1; % mark position with -2 if inside this query
            if mpchk > 0.6
                top3mp6 = top3mp6 + mpchk;  % more stringent criterion for m' only > 0.6 used.
                jktop36 = jktop36 + 1;
            end
        end
    end
    top6mpn = top6mpn/jk6;   % These average values of m' are without regard to size of m'
    if top3mpn > 0
        top3mpn = top3mpn/3;
    else
        top3mpn = top6mpn;
    end
    if jktop36 == 3
        top3mp6 = top3mp6/jktop36;
        top3d3 = top3d3 + 1;
    elseif jktop36 == 2 
        top3mp6 = top3mp6/jktop36;
        top3d2 = top3d2 + 1;
    elseif jktop36 == 1
        top3mp6 = top3mp6/jktop36;
        top3d1 = top3d1 + 1;
    elseif jktop36 == 0
        top3mp6 = 0;
        top3d0 = top3d0 + 1;
    end
    
    damp(:,gbnum) = [mploc mpmax 0 dpsum damage gbodam caxmis dpsum4 damage4 gbodam4 normdp normdp4 normp normp4 ...
        maxF1 maxF14 maxF1*Eratio maxF14*Eratio Eratio top3mpn top6mpn top3mp6]; % this is a summary matrix used for plotting
end  % of gbnum loop

fprintf(' %d\n ', gbnum);
% for i = 1:1:dRBdy   %  make smaller sorted version
%     sortmpr7(:,:,i) = mpr(1:7,1:7,i);
% end
plotname = {'      ','      ','      ','  dm''sum  ','  m*dm''sum  ','  gbo*m*dm''sum  ','  cax-mis  ','  dm''sum4  ',...
            '  m*dm''4 ','  gbo*dam4  ','  norm dm''  ','  norm dm''4  ','  norm m''  ','  norm m''4  ',...
            '  max(F1A,F1B)  ','  max(F14A,F14B)  ','  Emax(F1A,F1B)  ','  Emax(F14A,F14B)  ',' Eratio  ','top3mpn','top6mpn','top3mp6'}; 
%              4 = dm'sum   5 = m*dm'sum   6 = gbo*dam   7 = caxmis  8 = dm'sum4  
%              9 = damage4   10 = gbo*dam4   11 = norm dm'  12 = norm dm'4   13 = norm m'  14 = norm m'4 
%              15=max(F1A,F1B) 16=max(F14A,F14B) 17=Emax(F1A,F1B),
%              18=Emax(F14A,F14B), 19 = Eratio, 20 = top3mpn, 21 = top6mpn, 22 = top3mp6
% Definitions of parameters computed:  
%  dampar = dm' = adjustment of m' such that m' = 0.8 is worst case damage condition (i.e. dm' = 1 when m' = 0.8; When m' = 0.6 or 1, dm' = 0.8) 
%       This assumes that damage dm' would less if m' is higher because more slip transfer occurs with less debris left in the GB,
%       and that damage would be lower for smaller m' values because less slip transfer would happen
% 2 mpmax maximum value of m'
% 4 dpsum or dm'sum = sum of damage parameter dm' for all m' values computed, e.g. for m' values > 0.6 and slip systems having 
%       Schmid factors > tolL (lowest tolerance of Schmid factors, e.g. 0.2)
% 5 damage or m*dm'sum = Schmid factor times sum of dm' values computed, assuming that Schmid factor will scale the magnitude of slip transfer that occurs.
% 6 gbodam or gbo*m*dm''sum = damage parameter modified by apparent GB inclination (boundary normals || tensile axis see more damage).
% 7 caxmis or c-axis-misorientation is the angle between c-axes of the two grain orientations
% 8 dpsum4 or dm''sum4  Same as dpsum, but only for slip systems with m > tolH (e.g. 0.4)
% 9 damage4 or m*dm''sum4  Same as damage, but only for slip systems with m > tolH (e.g. 0.4)
% 10 gbodam4 or gbo*dam4  Same as gbodam, but only for slip systems with m > tolH (e.g. 0.4)
% 11 normdp   norm dm''   dpsum averaged (normalized by number of instances)
% 12 normdp4  norm dm''4  dpsum4 averaged (normalized by number of instances)
% 13 normp    norm m''    mpsum averaged (normalized by number of instances)
% 14 normp4   norm m''4   mpsum4 averaged (normalized by number of instances)
% 15 maxF1  the fip F1 is computed for all combinations for which m > tolL, and the largest value for grain A or B is chosen
% 16 maxF14 the fip F14 is computed for all combinations for which m > tolH, and the largest value for grain A or B is chosen
% 17 maxF1*Eratio  F1 times Eratio
% 18 maxF14*Eratio F14 times Eratio
% 19 Eratio  Ratio of Moduli such that it is < 1
% 20 top3mpn  Average of the m' values for the top 3 (or fewer) Schmid factor pairs in 5x5 upper left corner of the m' array
% 21 top6mpn  Average of the m' values for the top 6 (or fewer) Schmid factor pairs in 5x5 upper left corner of the m' array
% 22 top3mp6  Average of top 3 or fewer m' values > 0.6 in 5x5 upper left corner of the m' array


% % Choose which 4 histograms to plot and number of bins for each


bins = [25 25 25 25];  
pl = [4 5 6 7 8 9 10 19];        % Group A  or
pl = [13 14 11 12 15 16 17 18];  % Group B
figure, hold on 
for i = 1:1:4
    subplot(4,1,i), hist(damp(pl(i),:), bins(i)), ...
    title([num2str(min(damp(pl(i),:))), plotname{pl(i)}, num2str(bins(i)), ' bins, max = ', num2str(max(damp(pl(i),:)))]);
end
figure, hold on 
for i = 1:1:4
    subplot(4,1,i), hist(damp(pl(i+4),:), bins(i)), ...
    title([num2str(min(damp(pl(i+4),:))), plotname{pl(i+4)}, num2str(bins(i)), ' bins, max = ', num2str(max(damp(pl(i+4),:)))]);
end


%%  Choose what to plot along the grain boundary map from variables 4-19 noted above


plist = [ 20 21 22]; %13 14 11 12 15 16 17 18 19
for k = 1:1:length(plist) % 1:1:length(plist); %  4:4:4 %
    plnx = plist(k);
    if plnx == 13 || plnx == 14 
        bins = [.60 .64 .68 .72 .76 .80 .84 .88 .92 .96]; % for mp or mp4 or Eratio
    elseif plnx == 11 || plnx == 12
        bins = [.80 .82 .84 .86 .88 .90 .92 .94 .96 .98]; % for normp or normp4 
    elseif plnx == 7
        bins = [0 0.156 0.309 0.454 0.588 0.707 0.809 0.891 0.951 0.988 1.000];  % for c-axis plot    
    %   Corresponding degrees for bins above = [0 9 18 27 36 45 54 63 72 81 90];  % for c-axis plot
    else 
        bins = linspace(min(nonzeros(damp(plnx,:))), max(damp(plnx,:)), 11); % auto bin size for anything
    end
    
    binsdat = bins;   % temporary storage while setting up plot scale

    max1 = max(RBdy(:,9));         % max x value for grA
    max2 = max(RBdy(:,11));        % max x value for grB
    maxx = max([max1; max2]);
    may1 = max(RBdy(:,10));        % max y value for grA
    may2 = max(RBdy(:,12));        % max y value for grB
    maxy = max([may1; may2]);
    figure, hold on                % plot is based upon TV rastering, as given by TSL
    axis([0 maxx*2 -maxy .15*maxy ]); axis image ; %equal; 
    set(gcf,'Color',[.9,1,.9])      %  surrounding field is this color
    text(maxx*0.8,  0.08*maxy, plotname(plnx));
    text(maxx*-0.03, 0.18*maxy, num2str(bins,2));

    bins = [10 20 30 40 50 60 70 80 90 100];  % To set color key for boundaries
    for iGB = 1:1:100
        wid = 10;
        if fix(iGB/10)<iGB/10 % tick marks for color bar
            wid = 3;
        end
        widk = 1;    % thickness for width bar
        if iGB > 90
            widk=3;
        elseif iGB > 80
            widk=2;
        end
        vec = iGB;
        if vec<bins(1)                    %[.3 .3 .3] gray 0-bins(1)
            vRGB=[.3+.5*vec/bins(1)   .3-.3*vec/bins(1)  .3+.5*vec/bins(1)];               % gray -> purple                          
        elseif vec>=bins(1) && vec <bins(2);     %[.8 .0 .8] purple bins(1)-bins(2)
            vRGB=[.8-.8*(vec-bins(1))/(bins(2)-bins(1))  0.  .8-.1*(vec-bins(1))/(bins(2)-bins(1))];      % purple -> blue 
        elseif vec>=bins(2) && vec <bins(3);     %[.0 .0 .7] blue bins(2)-bins(3)
            vRGB=[0.  .9*(vec-bins(2))/(bins(3)-bins(2))  .7+.2*(vec-bins(2))/(bins(3)-bins(2))];         % blue -> turquoise    
        elseif vec>=bins(3) && vec <bins(4);     %[.0 .9 .9] turquoise bins(3)-bins(4)
            vRGB=[.1*(vec-bins(3))/(bins(4)-bins(3))  .9-.4*(vec-bins(3))/(bins(4)-bins(3))  .9-.8*(vec-bins(3))/(bins(4)-bins(3))]; % turquoise -> dk grn 
        elseif vec>=bins(4) && vec <bins(5);     %[.1 .5 .0] dk grn bins(4)-bins(5)
            vRGB=[.1+.1*(vec-bins(4))/(bins(5)-bins(4))  .5+.4*(vec-bins(4))/(bins(5)-bins(4))  .0];    % dk grn -> green 
        elseif vec>=bins(5) && vec <bins(6);     %[.2 .9 .0] green bins(5)-bins(6)
            vRGB=[.2+.7*(vec-bins(5))/(bins(6)-bins(5))  .9  0.];                          % green -> yellow                  
        elseif vec>=bins(6) && vec <bins(7);     %[.9 .9 .0] yellow bins(6)-bins(7)
            vRGB=[.9  .9-.3*(vec-bins(6))/(bins(7)-bins(6))  0.];                       % yellow -> orange                     
        elseif vec>=bins(7) && vec <bins(8);     %[.9 .6 .0] orange bins(7)-bins(8)                            
            vRGB=[.9+.1*(vec-bins(7))/(bins(8)-bins(7))  .6-.6*(vec-bins(7))/(bins(8)-bins(7))  0.]; % orange -> red
        elseif vec>=bins(8) && vec <bins(9);     %[1. .0 .0] red bins(8)-bins(9)                            
            vRGB=[1.  .7*(vec-bins(8))/(bins(9)-bins(8)) .7*(vec-bins(8))/(bins(9)-bins(8))]; % red -> pink
        elseif vec>=bins(9) && vec <bins(10);     %[1. .7 .7] pink bins(9)-bins(10)                            
            vRGB=[1.-.4*(vec-bins(9))/(bins(10)-bins(9))  .7-.7*(vec-bins(9))/(bins(10)-bins(9)) .7-.4*(vec-bins(9))/(bins(10)-bins(9))]; % red -> pink
        elseif vec > bins(10);
            vRGB = [.6  0. .4];                         %[.6  0. .3] mauve bins(9)
        end
    %    plot([1180;1180],maxx*[iGB-1; iGB],'Linewidth',(minF1m+iGB*dF1m-.95*minF1m)/(maxF1m-minF1m)*3,'color',[.1,.1,.1]);
    plot(maxx/150*[iGB-1; iGB],[0.1*maxy; 0.1*maxy],'Linewidth',widk,'color',[.1,.1,.1]);  % width bar plotted on top
    plot(maxx/150*[iGB-1; iGB],[0.05*maxy; 0.05*maxy],'Linewidth',wid,'color',vRGB);  % place color key x coordinate at suitable place...
    end

    bins = binsdat;  % replacing the necessary bins parameters for plotting along grain boundaries
    for iGB = 1:1:dRBdy(1,1)     % Plot parameter values on grain boundaries 
    %    wid = (RBdy(iGB,16)-.95*minF1m)/(maxF1m-minF1m)*4;
        wid = 2;
        if damp(plnx,iGB) > bins(9)
            wid = 4;
        elseif damp(plnx,iGB) > bins(8)
            wid = 3;
        end

        vec = damp(plnx,iGB);
        if vec<bins(1);                     %[.3 .3 .3] gray 0-10
            vRGB=[.3+.5*vec/bins(1)   .3-.3*vec/bins(1)  .3+.5*vec/bins(1)];               % gray -> purple                          
        elseif vec>=bins(1) && vec <bins(2);     %[.8 .0 .8] purple 10-20
            vRGB=[.8-.8*(vec-bins(1))/(bins(2)-bins(1))  0.  .8-.1*(vec-bins(1))/(bins(2)-bins(1))]; % purple -> blue 
        elseif vec>=bins(2) && vec <bins(3);     %[.0 .0 .7] blue 20-30
            vRGB=[0.  .9*(vec-bins(2))/(bins(3)-bins(2))  .7+.2*(vec-bins(2))/(bins(3)-bins(2))]; % blue -> turquoise    
        elseif vec>=bins(3) && vec <bins(4);     %[.0 .9 .9] turquoise 30-40
            vRGB=[.1*(vec-bins(3))/(bins(4)-bins(3))  .9-.4*(vec-bins(3))/(bins(4)-bins(3))  .9-.8*(vec-bins(3))/(bins(4)-bins(3))]; % turquoise -> dk grn 
        elseif vec>=bins(4) && vec <bins(5);     %[.1 .5 .0] dk grn 40-50
            vRGB=[.1+.1*(vec-bins(4))/(bins(5)-bins(4))  .5+.4*(vec-bins(4))/(bins(5)-bins(4))  0.];    % dk grn -> green 
        elseif vec>=bins(5) && vec <bins(6);     %[.2 .9 .0] green 50-60
            vRGB=[.2+.7*(vec-bins(5))/(bins(6)-bins(5))  .9  0.];                          % green -> yellow                  
        elseif vec>=bins(6) && vec <bins(7);     %[.9 .9 .0] yellow 60-70
            vRGB=[.9  .9-.3*(vec-bins(6))/(bins(7)-bins(6))  0.];                       % yellow -> orange                     
        elseif vec>=bins(7) && vec <bins(8);     %[.9 .6 .0] orange 70-80                            
            vRGB=[.9+.1*(vec-bins(7))/(bins(8)-bins(7))  .6-.6*(vec-bins(7))/(bins(8)-bins(7))  0.]; % orange -> red
        elseif vec>=bins(8) && vec <bins(9);     %[1. .0 .0] red 80-90                            
            vRGB=[1.  .7*(vec-bins(8))/(bins(9)-bins(8))  .7*(vec-bins(8))/(bins(9)-bins(8))]; % red -> pink
        elseif vec>=bins(9) && vec <bins(10);     %[1. .7 .7] pink bins(9)-bins(10)                            
            vRGB=[1.-.4*(vec-bins(9))/(bins(10)-bins(9))  .7-.7*(vec-bins(9))/(bins(10)-bins(9)) .7-.4*(vec-bins(9))/(bins(10)-bins(9))]; % red -> pink
        elseif vec > bins(10);
            vRGB = [.6  0. .3];                         %[.6  0. .4] mauve bins(9)
        end
        plot([RBdy(iGB,9);RBdy(iGB,11)],[-RBdy(iGB,10);-RBdy(iGB,12)],'Linewidth',wid,'color',vRGB); 
    end
end % k for plist


% %  plot grain numbers


for ng = 1:1:grmax 
	if grcen(ng,1) == 1 || grcen(ng,1) == 4 
        text(grcen(ng,2),-grcen(ng,3),mat2str(ng), 'FontWeight', 'bold');  %[ng grcen(ng,4:6)]
    elseif grcen(ng,1) == 2
        text(grcen(ng,2),-grcen(ng,3),int2str(ng),'color',[0 0 .5], 'FontWeight', 'bold');
    elseif grcen(ng,1) == 3
        text(grcen(ng,2),-grcen(ng,3),int2str(ng),'color',[0 .5 0], 'FontWeight', 'bold');
    end
end


% %  plot gb numbers


for gbnum = 1:1:dRBdy   % plot grain boundary numbers
    text((2*RBdy(gbnum,9)+RBdy(gbnum,11))/3,-(2*RBdy(gbnum,10)+RBdy(gbnum,12))/3,int2str(gbnum),...
        'color',[.8 0 0], 'FontWeight', 'light');
end

%%                  plot plane traces at grain center


legbnumt = 5;   %  This is the length of the plane traces in microns (units used for RBdy and IDgr)
% red for prism   % blue for basal   % yellow for pyr <c+a>   % green for pyr <a>
for ng = 1:1:grmax   
    if grcen(ng,1) == 1   % i.e. hexagonal
        tprm = 3; tbas = 3; tpya  = 2; tpyc = 2; ttw = 3;
        for k = 1:1:nslp(1)
            if sortmv(k,1,ng) > 33     % twins
                if ttw > 0 && sortmv(k,3,ng)> 0.4
                    vRGB = [.5 .5 0]; linw = ttw; ttw = ttw - 1;
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    plot(xs,ys,'Linewidth',linw,'color',vRGB);
                end
            elseif sortmv(k,1,ng) > 15 && sortmv(k,1,ng) < 34    % pyramidal c+a slip
                if tpyc > 0 && sortmv(k,3,ng)> 0.4
                    vRGB = [.8 .8 0]; linw = tpyc; tpyc = tpyc - 1;
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    plot(xs,ys,'Linewidth',linw,'color',vRGB);
                end
            elseif sortmv(k,1,ng) > 9 && sortmv(k,1,ng) < 16   % pyramidal a slip
                if tpya  > 0 && sortmv(k,3,ng)> 0.4
                    vRGB = [.2 .8 0]; linw = tpya; tpya = tpya  - 1;
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    plot(xs,ys,'Linewidth',linw,'color',vRGB);
                end
            elseif sortmv(k,1,ng) > 3 && sortmv(k,1,ng) < 10   % prism slip
                if tprm > 0 && sortmv(k,3)> 0.2
                    vRGB = [1 0 0]; linw = tprm; tprm = tprm - 1;
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    plot(xs,ys,'Linewidth',linw,'color',vRGB);
                end
            elseif sortmv(k,1,ng) < 4   % basal slip
                if tprm > 0 && sortmv(k,3,ng)> 0.2
                    vRGB = [0 0 1]; linw = tbas; tbas = tbas - 1;
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    plot(xs,ys,'Linewidth',linw,'color',vRGB);

                end
            end
        end
    elseif grcen(ng,1) == 2   % i.e. BCC   purple is (110), gold is (112) planes
        t110  = 3; t112 = 3; 
        for k = 1:1:nslp(2)
            if sortmv(k,1,ng) > 12     % 112 slip
                if t112 > 1 && sortmv(k,3,ng)> 0.3
                    vRGB = [.4 .2 .7];  % purple
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    if sortmv(k,3,ng)> 0.4
                        plot(xs,ys,'Linewidth',t112,'color',vRGB);
                    else
                        plot(xs,ys,':', 'Linewidth',t112-1,'color',vRGB);
                    end
                    t112 = t112 - 1;
                end  % for 112
            elseif sortmv(k,1,ng) < 13   % 110 slip
                if t110  > 1 && sortmv(k,3,ng)> 0.3
                    vRGB = [.7 .6 .1];   % gold
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    if sortmv(k,3,ng)> 0.4
                        plot(xs,ys,'Linewidth',t110,'color',vRGB);
                    else
                        plot(xs,ys,':', 'Linewidth',t110-1,'color',vRGB);
                    end
                    t110  = t110  - 1;
                end  % for 110
            end  % if 112 or 110
        end  % k
    elseif grcen(ng,1) == 3   % i.e. FCC   111 slip planes are brown
        t111  = 3;  
        for k = 1:1:nslp(3)
            if t111 > 1 && sortmv(k,3,ng)> 0.3
                vRGB = [0.3 0.3 0];  % brown
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                if sortmv(k,3,ng)> 0.4
                    plot(xs,ys,'Linewidth',t111,'color',vRGB);
                else
                    plot(xs,ys,':', 'Linewidth',t111-1,'color',vRGB);
                end
                t111 = t111 - 1;
            end  % for 111
        end  % k
    elseif grcen(ng,1) == 4   % i.e. BCT   purple is (110), gold is (112) planes
        tbe = 3; tbh = 3; t211 = 3; 
        for k = 1:1:nslp(4)
            if sortmv(k,1,ng) > 24     % 211 slip
                if t211 > 1 && sortmv(k,3,ng)> 0.3
                    vRGB = [.9 0 .8];  % magenta
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    if sortmv(k,3,ng)> 0.4
                        plot(xs,ys,'Linewidth',t211,'color',vRGB);
                    else
                        plot(xs,ys,':', 'Linewidth',t211-1,'color',vRGB);
                    end
                    t211 = t211 - 1;
                end  % for 211
            elseif sortmv(k,1,ng) < 13   % easier slip systems
                if tbe > 0 && sortmv(k,3,ng)> 0.4
                    vRGB = [.2 .8 .2];   % green
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    plot(xs,ys,'Linewidth',tbe,'color',vRGB);
                    tbe = tbe - 1;
                end
            else   % harder slip systems
                if tbh  > 0 && sortmv(k,3,ng)> 0.4
                    vRGB = [1 .6 .2];   % tan orange
                    if sortmv(k,11,ng) < 0
                        xs = [grcen(ng,2);grcen(ng,2)+legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)-legbnumt*sortmv(k,10,ng)];
                    else
                        xs = [grcen(ng,2);grcen(ng,2)-legbnumt*sortmv(k,11,ng)];
                        ys = [-grcen(ng,3);-grcen(ng,3)+legbnumt*sortmv(k,10,ng)];
                    end
                    plot(xs,ys,'Linewidth',tbh ,'color',vRGB);
                    tbh  = tbh  - 1;
                end
            end  % if 211 or other systems
        end  % k
    end  % if tree for crystal structure
end  % for ng (grain counter)


%%  Choose your favorite grain boundary ---------------------------------------
gbnum = GBNumChoice; 
whiteannotation = 1;  % make gb trace and axes white
ptpl = 0; % plots plane traces if = 1    Red dashed line in plot is perpendicular to the line connecting the centers of the two grains
prsxplt = zeros(8,7);    %               Black solid line in plot is the RC boundary segment orientation
prsyplt = zeros(8,7);     %  13 and 14 lead to misplaced plane traces.
prszplt = zeros(8,7); 
mp4gr = himp4(:,:,gbnum);
sorthimp4 = sortrows(mp4gr,-4);  % This sorts on basis of actual m' value

gbcen = [RBdy(gbnum,11)+RBdy(gbnum,9) RBdy(gbnum,12)+RBdy(gbnum,10)]/2; 
posv13c = [grcen(RBdy(gbnum,13),2) grcen(RBdy(gbnum,13),3)] - gbcen;  % find vector from center of GB to grain center 
posv14c = [grcen(RBdy(gbnum,14),2) grcen(RBdy(gbnum,14),3)] - gbcen;
v1314b = [(RBdy(gbnum,12)-RBdy(gbnum,10)) -(RBdy(gbnum,11)-RBdy(gbnum,9)) ]; % vector [dy,-dx] pointing perpendicular to GB
v1314a = v1314b;
if norm(posv13c)>norm(posv14c)
    if v1314b * posv14c' > 0   % use 14 to define negative vector for GB normal
        v1314b = -v1314b;       % This will make v1314c point to 13 when 14 is closer.
    end
else
    if v1314b * posv13c' < 0  % use 13 to define vector for GB normal
        v1314b = -v1314b;       % This will make v1314c point to 13 when 13 is closer.  
    end
end
if v1314b(1,1)<0  % determine if 13 is on the right or left
    Or13 = 'L g';   
    Or14 = 'R g';   
    gngbgn = num2str([RBdy(gbnum,13) gbnum RBdy(gbnum,14)]);
else
    Or13 = 'R g';
    Or14 = 'L g';   
    gngbgn = num2str([RBdy(gbnum,14) gbnum RBdy(gbnum,13)]);
end
v1314b = v1314b/norm(v1314b);  % find unit vector pointing perpendicular to GB in raster coords

ipl = -6; %  Strategy: Next, start isc loop for plotting slip systems
for imp = 1:1:size(sorthimp4,1)   % 1 imp is counter for m' values
    if sorthimp4(imp,1) ~= 0   %  evaluate only for recorded values (m'>.6)   
        if ipl-imp==-7 % six plots on a page
            figure
            ipl=ipl+6;
        end
        subplot(2,3,imp-ipl); hold on
        %figure;  hold on;
        axis square;
        axis([-3 3 -3 3]);     % plot representation of grain boundary
        if whiteannotation == 1
            set(gca ,'ycolor' ,'w'); set(gca ,'xcolor' ,'w');  % make axes white for ease in later arranging.
        else
            plot([0 1.5*cosd(RBdy(gbnum,8))], [0 1.5*sind(RBdy(gbnum,8))], '-k'); % plots gb from map from angle given in normal x-y space
        end
    %	text(0, 0, ['m''=', num2str(sorthimp4(imp,3),3)]);
        for igr = 13:1:14   % 1 i.e. first for grain in column 13, then 14 in Reconstructed Boundary file.
            if igr == 13
                issr = sorthimp4(imp,1);  % ss rank # in gr13   issr is slip system Schmid factor order # 
                grnum = RBdy(gbnum,13);
                del = v1314b;  
                cellcolor = [0 0 0]; %[.3 0 .5];                plot ([0 del(1)], -[0 del(2)])
            else     %  del is position vector from center of gb to 13 in raster coord system (ydown)
                grnum = RBdy(gbnum,14);
                issr = sorthimp4(imp,2);  % ss rank # in gr14  
                del = -v1314b;  
                cellcolor = [0 0 0]; %[0.5 0 0];                plot ([0 del(1)], -[0 del(2)])
            end
            ph = grcen(grnum,1);

                % (This orthorhombic (hexahedral cell) version was developed at a different time than the hex version, so details and approach differ)
            if ph > 1   %  plot the image of cubic unit cell,  slip vectors, planes, plane normals, and plane traces
                % sortmv contains vectors sorted to the bottom row (nslp+1) that define the unit cell
                x = [sortmv(mnslp+1,4:6,grnum)]; 
                y = [sortmv(mnslp+1,7:9,grnum)];
                z = [sortmv(mnslp+1,10:12,grnum)];
                ucell = [x;y;z]; % diagnostic
                % set up the 12 unit cell edge vectors    Based on X down and Y right
                prx = [0 0 0 x(1) x(2) x(3); ...
                    y(1) y(2) y(3) y(1)+x(1) y(2)+x(2) y(3)+x(3);...
                    z(1) z(2) z(3) z(1)+x(1) z(2)+x(2) z(3)+x(3);...
                    y(1)+z(1) y(2)+z(2) y(3)+z(3) y(1)+z(1)+x(1) y(2)+z(2)+x(2) y(3)+z(3)+x(3)];
                pry = [0 0 0 y(1) y(2) y(3); ...
                    x(1) x(2) x(3) y(1)+x(1) y(2)+x(2) y(3)+x(3);...
                    z(1) z(2) z(3) z(1)+y(1) z(2)+y(2) z(3)+y(3);...
                    x(1)+z(1) x(2)+z(2) x(3)+z(3) y(1)+z(1)+x(1) y(2)+z(2)+x(2) y(3)+z(3)+x(3)];
                prz = [0 0 0 z(1) z(2) z(3); ...
                    x(1) x(2) x(3) z(1)+x(1) z(2)+x(2) z(3)+x(3);...
                    y(1) y(2) y(3) z(1)+y(1) z(2)+y(2) z(3)+y(3);...
                    y(1)+x(1) y(2)+x(2) y(3)+x(3) y(1)+z(1)+x(1) y(2)+z(2)+x(2) y(3)+z(3)+x(3)];
                % figure out which are to be drawn (don't draw the lowest set in z direction)
                Sortprx = sortrows(prx,-3);
                Sortpry = sortrows(pry,-3);
                Sortprz = sortrows(prz,-3);
                % find range of x and y for plot scales
                prx1 = sortrows(prx,1);
                prx2 = sortrows(prx,4);
                pry1 = sortrows(pry,2);
                pry2 = sortrows(pry,5);
                minx = min(prx1(1,1), prx2(1,4));
                miny = min(pry1(1,2), pry2(1,5));
                maxx = max(prx1(4,1), prx2(4,4));
                maxy = max(pry1(4,2), pry2(4,5));
                cellcenter = [(minx+maxx)/2 (miny+maxy)/2];
                dx = del(2) - cellcenter(1);  % originally dx = del(1) - cellcenter(1); 
                dy = del(1) - cellcenter(2);  % originally dy = del(2) - cellcenter(2);
                [dx dy]; % diagnostic               

                p1 = sortmv(issr,13:15,grnum);       % points on the slip plane
                p2 = sortmv(issr,16:18,grnum);
                p3 = sortmv(issr,19:21,grnum);
                p4 = sortmv(issr,22:24,grnum);
                p5 = sortmv(issr,25:27,grnum);
                p6 = sortmv(issr,28:30,grnum);
                spx = [p1(1) p2(1) p3(1) p4(1) p5(1) p6(1)]+dx;
                spy = [p1(2) p2(2) p3(2) p4(2) p5(2) p6(2)]+dy;
                ssn = sortmv(issr,1,grnum);           % slip system number
                Sf = sortmv(issr,2,grnum);           % Schmid factor
    %            n = [p1 p3]; % p1+sortmv(issr,4:6,grnum)];  % plane normal
                b = [p1 p3]; % p1+sortmv(issr,7:9,grnum)];  % Burgers vector
                pt = sortmv(issr,10:12,grnum);       % plane trace

            % These plots will match TSL with X down !!!!

                plot([0+dy x(2)+dy],-[0+dx x(1)+dx], ':', 'Linewidth',3,'Color',[1 0 0]);% plot x = red
                plot([0+dy y(2)+dy],-[0+dx y(1)+dx], ':', 'Linewidth',3,'Color',[.6 .7 0]);% plot y = green-gold
                plot([0+dy z(2)+dy],-[0+dx z(1)+dx], ':', 'Linewidth',3,'Color',[0 0 1]);% plot z = blue

    % % 	   %text(min_x, min_y-1, 'Eulers = ',num2str(grcen(gpl(igr-1),4)), num2str(grcen(gpl(igr-1),5)), num2str(grcen(gpl(igr-1),6)));

                if sortmv(issr,6,grnum)>0
                     fill(spy,-spx, [.8 .8 .65])  % slip plane filled warm gray if top side
    %                 plot([n(2) n(5)]+dy, -([n(1) n(4)]+dx),'Linewidth',3,'Color',[.8 .8 .65]);
                else                
                     fill(spy,-spx, [.65 .65 .7])  % slip plane filled cool gray if bottom side
    %                 plot([n(2) n(5)]+dy, -([n(1) n(4)]+dx),'Linewidth',3,'Color',[.65 .65 .7]);
                end
                if sortmv(issr,6,grnum)>0
                    Bvcolor = [0 .7 .7];
                else
                    Bvcolor = [0 1 1];
                end
                Sfs = sign(Sf); % diagnostic
                if Sf > 0    % plot Burgers vector direction using first and third points in slip plane
                    plot(b(2)+dy,-(b(1)+dx),'.','MarkerSize', 24, 'Color', Bvcolor)
                    plot([b(2) b(5)]+dy,-([b(1) b(4)]+dx),'Linewidth',3,'Color',Bvcolor)
                else         % plot Burgers vector in opposite direction
                    plot(b(5)+dy,-(b(4)+dx),'.','MarkerSize', 24, 'Color', Bvcolor)
                    plot([b(5) b(2)]+dy,-([b(4) b(1)]+dx),'Linewidth',3,'Color',Bvcolor)
                end

                for i = 1:1:3 % plot prism
                    plot([Sortprx(i,2)+dy Sortprx(i,5)+dy],-[Sortprx(i,1)+dx Sortprx(i,4)+dx], 'Linewidth' ,1.5, 'Color', cellcolor);
                    plot([Sortpry(i,2)+dy Sortpry(i,5)+dy],-[Sortpry(i,1)+dx Sortpry(i,4)+dx], 'Linewidth' ,1.5, 'Color', cellcolor);
                    plot([Sortprz(i,2)+dy Sortprz(i,5)+dy],-[Sortprz(i,1)+dx Sortprz(i,4)+dx], 'Linewidth', 1.5, 'Color', cellcolor);
                end
                if Sortprx(4,1) ~= 0  %  overwrite x-y-z axes unless they are underneath the plane
                    plot([0+dy x(2)+dy],-[0+dx x(1)+dx], ':', 'Linewidth',3,'Color',[1 0 0]);% plot x = red
                end
                if Sortpry(4,1) ~= 0 
                    plot([0+dy y(2)+dy],-[0+dx y(1)+dx], ':', 'Linewidth',3,'Color',[.6 .7 0]);% plot y = gold
                end
                if Sortprz(4,1) ~= 0 
                    plot([0+dy z(2)+dy],-[0+dx z(1)+dx], ':', 'Linewidth',3,'Color',[0 0 1]);% plot z = blue
                end
                if ptpl == 1
                    if ph == 2
                        if ssn<25 && ssn>12   % {112 slip directions}  purple
                            plot([-pt(2)+dy pt(2)+dy],-[-pt(1)+dx pt(1)+dx],'--','Linewidth',3,'Color',[.4 .2 .7]) 
                        else          %  {110 slip directions}  gold 
                            plot([-pt(2)+dy pt(2)+dy],-[-pt(1)+dx pt(1)+dx],'--','Linewidth',3,'Color',[.7 .6 .1]) 
                        end                %---->  NOTE that Schmid factor vector is plotted in correct direction, 
                    elseif ph == 3
                        plot([-pt(2)+dy pt(2)+dy], -[-pt(1)+dx pt(1)+dx],'--','Linewidth',2,'Color',[.3 .3 0])
                    elseif ph == 4
                        if ssn>24  % {211} plane traces      magenta
                            plot([-pt(2)+dy pt(2)+dy], -[-pt(1)+dx pt(1)+dx],'--','Linewidth',2,'Color',[.9 0 .8]) 
                        elseif ssn<13   % {easier slip planes} green
                            plot([-pt(2)+dy pt(2)+dy], -[-pt(1)+dx pt(1)+dx],'--','Linewidth',3,'Color',[.2 .8 .2]) 
                        else          %  {medium slip systems}  tan orange
                            plot([-pt(2)+dy pt(2)+dy], -[-pt(1)+dx pt(1)+dx],'--','Linewidth',2,'Color',[1 .6 .2]) 
                        end
                    end
                end
                if ph == 2
                    sliplane = ssbcc(1,:,ssn);  Burgvec = ssbcc(2,:,ssn);
                elseif ph == 3
                    sliplane = ssfcc(1,:,ssn);  Burgvec = ssfcc(2,:,ssn);
                elseif ph == 4
                    sliplane = ssbct(1,:,ssn);  Burgvec = ssbct(2,:,ssn);
                end
                plot(0+dy,-(0+dx),'ko');     %---->  and labeled with corretly signed b vector as plotted.
                if igr == 13
                    line1 = [Or13 num2str(grnum) ' m' num2str(issr) ' = ' num2str(Sfs*Sf, 2) '  ss' num2str(ssn)...
                        ' n' mat2str(sliplane)  mat2str(Sfs*Burgvec) 'b'];
                end
                if igr == 14
                    title({[line1] [Or14 num2str(grnum) ' m' num2str(issr) ' = ' num2str(Sfs*Sf, 2) '  ss' num2str(ssn)...
                    ' n' mat2str(sliplane)  mat2str(Sfs*Burgvec) 'b'] [gngbgn, ' m'' = ', num2str(sorthimp4(imp,4),3)] } );
                end


            else %----- i.e. for (grnum,1) == 1  for hexagonal  -------------
    %  Plot the image of hexagonal unit cell, slip vectors, planes, plane normals, and plane traces
    %  Strategy:  First extract useful vectors to draw the hexagonal prisms from slip system information 
    %     positions in mvs  p1:13-15  p2:16-18  p3:19-21  p4:22-24  p5:25-27  p6:28-30 
    %     positions in pln   p1:4-6    p2:7-9   p3:10-12  p4:13-15  p5:16-18  p6:18-21
                for isc = 1:1:nslphex 
                    if sortmv(isc,1,grnum) == 1;                % locate basal planes using SS1
                        pln(1,4:21) = sortmv(isc,13:30,grnum);  % bottom basal plane
                        pln(2,4:21) = sortmv(isc,13:30,grnum);  % top basal plane
                        rotc = sortmv(isc,4:6,grnum)*c_a_hex;       % basal plane normal * c/a
                        for j = 4:3:19
                            pln(2,j:j+2) = pln(1,j:j+2) + rotc;    % move top plane up by a unit of c
                        end
                        a1 = sortmv(isc,7:9,grnum);   %  locate a1 using SS1
                    elseif sortmv(isc,1,grnum) == 2
                        a2 = sortmv(isc,7:9,grnum);   %  locate a2 using SS2
                    elseif sortmv(isc,1,grnum) == 3
                        a3 = sortmv(isc,7:9,grnum);   %  locate a3 using SS3
                    end
                end
                for isc = 1:1:nslphex   
                    if sortmv(isc,1,grnum) == 4      %  locate two prism planes on opposite sides using SS4
                        pln(3,4:21) = sortmv(isc,13:30,grnum);
                        for j = 13:3:28
                            pln(4,j-9:j-7) = sortmv(isc,j:j+2,grnum) + a2 - a3;
                        end
                    elseif sortmv(isc,1,grnum) == 5;    %  locate two prism planes on opposite sides using SS5
                        pln(5,4:21) = sortmv(isc,13:30,grnum);
                        for j = 13:3:28
                            pln(6,j-9:j-7) = sortmv(isc,j:j+2,grnum) + a3 - a1;
                        end
                    elseif sortmv(isc,1,grnum) == 6;    %  locate two prism planes on opposite sides using SS6
                        pln(7,4:21) = sortmv(isc,13:30,grnum);
                        for j = 13:3:28
                            pln(8,j-9:j-7) = sortmv(isc,j:j+2,grnum) + a1 - a2;
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

                sp1 = sortmv(issr,13:15,grnum);       % identify plotted points on the slip plane
                sp2 = sortmv(issr,16:18,grnum);
                sp3 = sortmv(issr,19:21,grnum);
                sp4 = sortmv(issr,22:24,grnum);
                sp5 = sortmv(issr,25:27,grnum);
                sp6 = sortmv(issr,28:30,grnum);
                spx = [sp1(1) sp2(1) sp3(1) sp4(1) sp5(1) sp6(1) sp1(1)];
                spy = [sp1(2) sp2(2) sp3(2) sp4(2) sp5(2) sp6(2) sp1(2)];
                ssn = sortmv(issr,1,grnum);           % slip system number
                Sf = sortmv(issr,2,grnum);           % Schmid factor
                n = [0 0 0 sortmv(issr,4:6,grnum)];  % plane normal
                b = [sp1 sp4]; % p1+sortmv(issr,7:9,grnum)];  % Burgers vector
                pt = sortmv(issr,10:12,grnum);       % plane trace
    %             minx = min(minx, n(4));
    %             miny = min(miny, n(5));
    %             maxx = max(maxx, n(4));   % find appropriate range of x and y for plot to include plane normal
    %             maxy = max(maxy, n(5));
                midx = (minx+maxx)/2;
                midy = (miny+maxy)/2;
                cellcenter = [midx midy];
                dx =(1.6*del(2) - cellcenter(1));  % del is raster, cellcenter is TSL coords 
                dy =(1.6*del(1) - cellcenter(2));  % so not dy = del(2) - cellcenter(2);
                [dx dy]; % diagnostic               

            % These plots will match TSL with X down !!!!   Plotting starts
                if grcen(grnum,5) < 90 % if PHI < 90, then make the 3 coordinate axes visible below slip planes
                    plot([0 a1(2)]+dy,-([0 a1(1)]+dx), ':', 'Linewidth',3,'Color',[1 0 .2]);% plot x = red
                    plot([0 a2(2)]+dy,-([0 a2(1)]+dx), ':', 'Linewidth',3,'Color',[.6 .8 0]);% plot y = green-gold
                    plot([0 a3(2)]+dy,-([0 a3(1)]+dx), ':', 'Linewidth',3,'Color',[0 0 1]);% plot z = blue
                end

                if sortmv(issr,6,grnum)>0   % is k component of slip plane normal positive or negative?
                    fill(spy+dy,-(spx+dx), [.8 .8 .65])  % slip plane filled warm gray
            %        plot([n(2) n(5)], -[n(1) n(4)],'Linewidth',3,'Color',[.8 .8 .65]);
                else                % slip plane filled cool gray if normal has neg z component
                    fill(spy+dy,-(spx+dx), [.65 .65 .7])  
            %        plot([n(2) n(5)], -[n(1) n(4)],'Linewidth',3,'Color',[.65 .65 .7]);
                end
                if sortmv(issr,6,grnum)>0
                    Bvcolor = [0 .7 .7];
                    if ssn >= iC1
                        Bvcolor = [.1 .6 0];
                    end
                    if ssn >= iT1 && ssn <= fT2
                        Bvcolor = [1 .6 0];
                    end
                else
                    Bvcolor = [0 1 1];
                    if ssn >= iC1
                        Bvcolor = [.3 .9 0];
                    end
                    if ssn >= iT1 && ssn <= fT2
                        Bvcolor = [1 .8 0];
                    end
                end
                Sfs = 1;
                if ssn < iT1
                    Sfs = sign(Sf);
                end
                if Sf > 0    % plot Burgers vector direction
                    if ssn >= iT1    % this is for twins - the Burgers vector length is shown to be 1/2 of the usual length in the unit cell 
                        plot(b(2)+dy,-(b(1)+dx),'.','MarkerSize', 24, 'Color', Bvcolor)
                        plot([b(2) (b(2)+b(5))/2]+dy,-([b(1) (b(1)+b(4))/2]+dx),'Linewidth',4,'Color',Bvcolor)
                    else
                        plot(b(2)+dy,-(b(1)+dx),'.','MarkerSize', 24, 'Color', Bvcolor)
                        plot([b(2) b(5)]+dy,-([b(1) b(4)]+dx),'Linewidth',4,'Color',Bvcolor)
                    end
                else         % plot Burgers vector in opposite direction
                    if ssn >= iT1    % this is for twins - the Burgers vector length is shown to be 1/2 of the usual length in the unit cell 
                        plot(b(2)+dy,-(b(1)+dx),'.','MarkerSize', 24, 'Color', Bvcolor)
                        plot([b(2) (2*b(2)+b(5))/3]+dy,-([b(1) (2*b(1)+b(4))/3]+dx),'Linewidth',4,'Color',Bvcolor)
                    else
                        plot(b(5)+dy,-(b(4)+dx),'.','MarkerSize', 24, 'Color', Bvcolor)
                        plot([b(5) b(2)]+dy,-([b(4) b(1)]+dx),'Linewidth',4,'Color',Bvcolor)
                    end
                end

                for j = 1:1:4 % plot the 4 top most surface prisms of the hex cell that have the highest z elevation
                    plot(prsyplt(j,:)+dy,-(prsxplt(j,:)+dx), 'Linewidth', 2, 'Color', cellcolor);
                end
                if grcen(grnum,5) > 90 % if PHI < 90, make the 3 coordinate axes visible above slip planes
                    plot([0 a1(2)]+dy,-([0 a1(1)]+dx), ':', 'Linewidth',3,'Color',[1 0 .2]);% plot x = red
                    plot([0 a2(2)]+dy,-([0 a2(1)]+dx), ':', 'Linewidth',3,'Color',[.6 .8 0]);% plot y = green-gold
                    plot([0 a3(2)]+dy,-([0 a3(1)]+dx), ':', 'Linewidth',3,'Color',[0 0 1]);% plot z = blue
                end
                if ptpl == 1
                    if ssn >= iC1  
                        ptrcolor = [.2 .8 0]; % compression twin plane traces   green
                    elseif ssn >= iT1 && ssn <= fT2  
                        ptrcolor = [1 .6 0]; % extension twin plane traces   orange
                    elseif ssn > ipyrc && ssn <= f2pyrc   
                        ptrcolor = [.95 .85 0]; % <c+a> plane traces   green-gold
                    elseif ssn <= fpyra && ssn >= ipyra   
                        ptrcolor = [0 .9 .5]; % pyr <a>  green-blue
                    elseif ssn <= f2prs && ssn >= iprs    
                        ptrcolor = [1 .2 0]; % prism <a,aa>  red
                    else                             
                        ptrcolor = [0 0 1]; %  basal <a>  blue
                    end                %---->  NOTE that Schmid factor vector is plotted in correct direction, 
                    plot([-pt(2) pt(2)]+dy, -([-pt(1) pt(1)]+dx),'--','Linewidth',3,'Color',ptrcolor)
                end
                if igr == 13
                    line1 = [Or13 num2str(grnum) ' m' num2str(issr) ' = ' num2str(Sfs*Sf, 2) '  ss' num2str(ssn)...
                        ' n' mat2str(sshex(1,:,ssn))  mat2str(Sfs*sshex(2,:,ssn)) 'b'];
                end
                if igr == 14
                    title({[line1] [Or14 num2str(grnum) ' m' num2str(issr) ' = ' num2str(Sfs*Sf, 2) '  ss' num2str(ssn)...
                    ' n' mat2str(sshex(1,:,ssn))  mat2str(Sfs*sshex(2,:,ssn)) 'b'] [gngbgn, ' m'' = ', num2str(sorthimp4(imp,4),3)] } );
%                      text(-2.8, -2.8, ['Or-' num2str(gpl(igr-1)), mat2str([grcen(gpl(igr-1),4:6)]), ' Eulers ', mat2str(grcen(gpl(igr-1),4:6)), ...
%                      'Or-', num2str(grnum)] );
%                      text(minx, miny, ['Eulers = ',mat2str(grcen(gpl(igr-1),4:6))]);
                end
            end  % graincen for cubic or hex  if (has middle else)
        end   % igr - grains 1 and 2
    else
        fprintf('No slip system has Schmid factors > %4.2f for one of the two grains at GB %d %d\n', schmid_tolH, gbnum, imp)
    end  %  if statement for continuing the loop for non-zero himp4 values
end   % imp  m' loop









end