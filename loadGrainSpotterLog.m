function log = loadGrainSpotterLog(fileName,fPrefix)

% Load in a Grainspotter log file (based on parseGrainSpotterLog)
%
%   log = loadGrainSpotterLog(fileName) reads the Grainspotter log file
%   with the name fileName and returns the information in an array of
%   structures with fields:
%       nExpGvec = Number of expected G vectors
%       nMeasGvec = Number of measured G vectors
%       nMeasOnce = Number of G vectors measured once
%       nMeasMore = Number of G vectors measured more than once
%       meanIA = Average internal angle between prediced and measured
%       U = 3x3 Orientation matrix
%       gvec = G vector table
%       hkl = 3 hkl values
%
%   The columns of the log file are:
%     n gvector_id peak_id  h k l  h_pred k_pred l_pred  dh dk dl
%     tth_meas tth_pred dtth  omega_meas omega_pred domega
%     eta_meas  eta_pred deta  IA
%
%   Example:
%     log = parseGrainSpotterLog('simul.log');

fid=fopen(fullfile(fPrefix,fileName),'r');
%fid = fopen(fileName, 'r');

if(fid == -1)
    beep;
    error('Cannot open file:\n  %s\n', fileName);
end
% Get number of grains
line = fgetl(fid);
parms = sscanf(line, 'Found %i');
nGrains = parms;
textscan(fid, '%*[^\n]', 17); % Skip 22 lines
% Loop over found grains
log(nGrains) = struct('nExpGvec',[],'nMeasGvec',[],'nMeasOnce',[], ...
  'nMeasMore',[],'meanIA',[],'U',[],'UBI',[],'r',[],'euler',[],'q',[],'refl',[]);

for i = 1:nGrains
    fgetl(fid);  % Grain nnn, nPairs (Skip)
    
    parms = fscanf(fid, '%f', [1, 4]);              % #expected gvectors #measured gvectors #measured once #measured more than once
    
    log(i).nExpGvec = parms(1);
    log(i).nMeasGvec = parms(2);
    log(i).nMeasOnce = parms(3);
    log(i).nMeasMore = parms(4);   
    
    pos   = fscanf(fid, '%f', [1, 5]);              % mean_IA position_x position_y position_z pos_chisq
    log(i).meanIA = pos(1);
    
    U     = fscanf(fid, '%f', [3, 3]);              % U matrix
    log(i).U = U';
    
    Ub    = fscanf(fid, '%f', [3, 3]);               % UBI matrix
    log(i).UBI = Ub';
    
    r     = fscanf(fid, '%f', [1, 3]);              % r1 r2 r3
    log(i).r = r;
    
    euler = fscanf(fid, '%f', [1, 3]);              % phi1 phi phi2
    log(i).euler = euler;
    
    q     = fscanf(fid, '%f', [1, 4]);              % q0 qx qy qz
    log(i).q = q;
    
    C = textscan(fid, '%f', 22*log(i).nMeasGvec);
    % gvec is the whole 22 x nMeasGvec array of columns
    log(i).refl = reshape(C{1},22,log(i).nMeasGvec)';
    
    textscan(fid, '%*[^\n]', 1); % Skip 1 lines
end

fclose(fid);