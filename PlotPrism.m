function PlotPrism(Euler)

%%  plot the image of unit cell, slip vectors, planes, plane normals, and plane traces
iAng =1; %This value should be consistent with the number of the grain that is to be plotted!   
ptpl = 0;  % 0 = don't plot plane traces
%  Strategy:  First extract useful vectors from slip system information to draw the hexagonal prisms
%  positions in sortmv  p1:13-15  p2:16-18  p3:19-21  p4:22-24  p5:25-27  p6:28-30 
%  positions in pln    p1:4-6   p2:7-9   p3:10-12  p4:13-15  p5:16-18  p6:18-21
STens=[0,0,1,0,0,0];
nslphex=57;
c_a=1.587;
[sortmv,~]=HexEval(Euler,STens);
pln=zeros(8,21);
for isc = 1:1:nslphex   
    if sortmv(isc,1,iAng) == 1              % locate the two basal planes
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
    if sortmv(isc,1,iAng) == 4      %  locate two prism planes on opposite sides using SS4
    	pln(3,4:21) = sortmv(isc,13:30,iAng);
        for j = 13:3:28
            pln(4,j-9:j-7) = sortmv(isc,j:j+2,iAng) + a2 - a3;
        end
    elseif sortmv(isc,1,iAng) == 5   %  locate two prism planes on opposite sides using SS5
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
for isc = 1:1:1  % change nslphex to smaller number to plot fewer
    if ipl-isc==-9 % eight plots on a page
        figure
        ipl=ipl+8;
    end
    subplot(1,1,isc-ipl)
    hold on
    sp1 = sortmv(isc,13:15,iAng);       % identify plotted points on the slip plane
    sp2 = sortmv(isc,16:18,iAng);
    sp3 = sortmv(isc,19:21,iAng);
    sp4 = sortmv(isc,22:24,iAng);
    sp5 = sortmv(isc,25:27,iAng);
    sp6 = sortmv(isc,28:30,iAng);
    spx = [sp1(1) sp2(1) sp3(1) sp4(1) sp5(1) sp6(1) sp1(1)];
    spy = [sp1(2) sp2(2) sp3(2) sp4(2) sp5(2) sp6(2) sp1(2)];
    ssn = sortmv(isc,1,iAng);           % slip system number
    Sf = sortmv(isc,2,iAng);           % Schmid factor
    n = [0 0 0 sortmv(isc,4:6,iAng)];  % plane normal
    b = [sp1 sp4]; % p1+sortmv(isc,7:9,iAng)];  % Burgers vector
    pt = sortmv(isc,10:12,iAng);       % plane trace
    minx = min(minx, n(4));
    miny = min(miny, n(5));
    maxx = max(maxx, n(4));   % find appropriate range of x and y for plot
    maxy = max(maxy, n(5));
    midx = (minx+maxx)/2;
    midy = (miny+maxy)/2;
    del = 1.4;

% These plots will match TSL with X down !!!!   Plotting starts...
    axis square
    set(gca ,'ycolor' ,'w'); set(gca ,'xcolor' ,'w');  % make axes white for ease in later arranging.
    axis([midy-del midy+del -midx-del -midx+del])

    

    for j = 1:1:4 % plot the 4 top most surface prisms of the hex cell that have the highest z elevation
        plot(prsyplt(j,:),-prsxplt(j,:), 'Linewidth',2,'Color',[.0 .0 .0]);
    end

end





end