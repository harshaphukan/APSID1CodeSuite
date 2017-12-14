function [g0Cry,g_idl] = g0hcp( a0, c0)
%
% create hkl's
% pedestrian hard-coding approach. hkls taken from EQUIV using space group
% P63/mmc

hkil_{1}  = { 1  0 -1  0; 0  1 -1  0; 1 -1  0  0};
hkil_{2}  = { 0  0  0  2};
hkil_{3}  = { 1  0 -1  1; 0  1 -1  1; 1 -1  0  1; 1  0 -1 -1; 0  1 -1 -1; 1 -1  0 -1};
hkil_{4}  = { 1  0 -1  2; 0  1 -1  2; 1 -1  0  2; 1  0 -1 -2; 0  1 -1 -2; 1 -1  0 -2};
hkil_{5}  = { 2 -1 -1  0;-1  2 -1  0;-1 -1  2  0};
hkil_{6}  = { 1  0 -1  3; 0  1 -1  3; 1 -1  0  3; 1  0 -1 -3; 0  1 -1 -3; 1 -1  0 -3};
hkil_{7}  = { 2  0 -2  0; 0  2 -2  0; 2 -2  0  0};
hkil_{8}  = { 2 -1 -1  2; 2 -1 -1 -2;-1  2 -1  2;-1  2 -1 -2;-1 -1  2  2;-1 -1  2 -2};
hkil_{9}  = { 2  0 -2  1; 0  2 -2  1; 2 -2  0  1; 2  0 -2 -1; 0  2 -2 -1; 2 -2  0 -1};
hkil_{10} = { 0  0  0  4};
hkil_{11} = { 2  0 -2  2; 0  2 -2  2; 2 -2  0  2; 2  0 -2 -2; 0  2 -2 -2; 2 -2  0 -2};
hkil_{12} = { 1  0 -1  4; 0  1 -1  4; 1 -1  0  4; 1  0 -1 -4; 0  1 -1 -4; 1 -1  0 -4};
hkil_{13} = { 2  0 -2  3; 0  2 -2  3; 2 -2  0  3; 2  0 -2 -3; 0  2 -2 -3; 2 -2  0 -3};
hkil_{14} = { 2  1 -3  0; 2 -3  1  0;-3  2  1  0; 1  2 -3  0; 1 -3  2  0;-3  1  2  0};
hkil_{15} = { 2  1 -3  1; 2 -3  1  1;-3  2  1  1; 1  2 -3  1; 1 -3  2  1;-3  1  2  1; 2  1 -3 -1; 2 -3  1 -1;-3  2  1 -1; 1  2 -3 -1; 1 -3  2 -1;-3  1  2 -1 };
hkil_{16} = { 2 -1 -1  4; 2 -1 -1 -4;-1  2 -1  4;-1  2 -1 -4;-1 -1  2  4;-1 -1  2 -4};
% hkil_{17} = { 1  0 -1  5; 0  1 -1  5; 1 -1  0  5; 1  0 -1 -5; 0  1 -1 -5; 1 -1  0 -5};
% hkil_{18} = { 2  1 -3  2; 2 -3  1  2;-3  2  1  2; 1  2 -3  2; 1 -3  2  2;-3  1  2  2; 2  1 -3 -2; 2 -3  1 -2;-3  2  1 -2; 1  2 -3 -2; 1 -3  2 -2;-3  1  2 -2 };
% hkil_{19} = { 2  0 -2  4; 0  2 -2  4; 2 -2  0  4; 2  0 -2 -4; 0  2 -2 -4; 2 -2  0 -4};
% hkil_{20} = { 3  0 -3  0; 0  3 -3  0; 3 -3  0  0};

hkil = [];
  for i = 1:length(hkil_)
    hkil = [hkil; cell2mat(hkil_{i}); -cell2mat(hkil_{i})];
  end
hkl = hkil(:,[1 2 4])'; % 3 notation

% A matrix

% The lattice parameters:
% Joel xc || a, yc in (a, b) plane, zc = cross(xc, yc)
al = 90; be = 90; ga = 120; 
% from calling statement a0 = 2.937; b0 = a0; c0 = 4.678; % simplify
b0 = a0;

as = acosd(( cosd(be)*cosd(ga)-cosd(al) )/ sind(be) / sind(ga));
A0 = [a0 b0*cosd(ga)  c0*cosd(be); ...
      0  b0*sind(ga) -c0*sind(be)*cosd(as); ...
      0  0            c0*sind(be)*sind(as)];
  
% lam = 0.18972; % from call to function

g0Cry = A0'\hkl; % unstrained reciprocal lattice vectors

g_idl = zeros(length(hkil_),1);

lastlen = 0;
idx = 0;
for i=1:length(g0Cry)
    gl = norm(g0Cry(:,i));
     if ( abs(gl-lastlen) > 0.0001 )
         idx = idx + 1;
         g_idl(idx) = gl;
         lastlen = gl;
     end
         
end

g_idl = sort(g_idl);
