function [ gve ] = loadGVE( fileName )
%loadGVE load g vectors from .gve file

fid = fopen(fileName, 'r');
if(fid == -1)
    beep;
    error('Cannot open file:\n  %s\n', fileName);
end

flag = 0;

% Read through, until the listing of g-vectors
while ~flag

    line = fgetl(fid);
    
    flag = isequal(line,'#  gx  gy  gz  xc  yc  ds  eta  omega  spot3d_id  xl  yl  zl');
  
end

gve = fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g', [12 inf]);

gve = gve';

fclose(fid);

end

