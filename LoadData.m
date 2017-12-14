%% Function to load data from .log and .gve files and Output an Array with Grain number, Center-of-Mass-Coordinates, StressTensor and Euler Angles
% Note that the Euler angles are corrected by adding 30-degrees to phi2 in
% the Result
function [GrainArray]=LoadData()
%Routine to extract strain/stress tensors and positional data of grains
%from .gve and .log file. Developed by Beaudoin and Wang: Uses least
%squaresapproach proposed by Marguiles et al. [2002]

fPrefixGVE=input('Type in path for .gve file:','s');
fPrefixLOG=input('Type in path for .log file:','s');

file_nameLog=input('Enter log filename, e.g xyz.log:','s');
log_data=loadGrainSpotterLog(file_nameLog,fPrefixLOG);


file_nameGVE=input('Enter GVE filename, e.g xyz.gve:','s');

fPostfix = '';
[name,ext]=fileparts(file_nameLog);
logStem=strtok(ext,'_');
fLogStem=strcat(logStem,'_','00');


fGVEstem=fLogStem;

[nameG,extG]=fileparts(file_nameGVE);

[a,fPostfixGVE]=strtok(extG,'t');
fPostfixGVE=strcat('_',fPostfixGVE);

[b,fPostfixLOG]=strtok(ext,'t');
fPostfixLOG=strcat('_',fPostfixLOG);



N_grain=length(log_data); 

%% This part extracts the Image Number from the .log filename
if length(ext)==22
    Image_string=strcat(ext(9),ext(10),ext(11));
elseif length(ext)==21
    Image_string=strcat(ext(8),ext(9),ext(10));
else
    print('Please rename your .log and .gve files to fit the compatible format.\n')
end
if Image_string(1)=='0'
    ImageString=strcat(Image_string(2),Image_string(3));
else
    ImageString=Image_string;
end

Image_number=str2num(ImageString);

scanNo = Image_number*ones(N_grain,1);

grainNo = 1:N_grain;
 
strain=zeros(3,3);


%[g0Cry, g_idl] = g0hcp(2.937, 4.678); % Armand
[g0Cry, g_idl] = g0hcp(2.95, 4.683); % New parameter

  T1twin=[0.586, 0.338, 0.737, -0.638, -0.368, 0.676; 0, 0.676, 0.737, 0, -0.737, 0.676; -0.586, 0.338, 0.737, 0.638, -0.368, 0.676;-0.586, -0.338, 0.737, 0.638, 0.368, 0.676; 0 -0.676 0.737 0 0.737 0.676; 0.586 -0.338 0.737 -0.638 0.368 0.676]; % twin plane and twin direction in the crystal coordinate system for the six T1 twin varients  
    for tn=1:6
    mT1twin(:,:,tn)=(kron(T1twin(tn, 1:3)',T1twin(tn, 4:6))+kron(T1twin(tn, 4:6)',T1twin(tn, 1:3)))/2; % dyadic product for later calculating resolved shear stress
    end
  
g_ideal = g_idl;
EulerAng=zeros(length(scanNo),3);
GrainArray=zeros(length(scanNo),12);

for scan = 1 : length(scanNo)
    
    filelog = parseGrainSpotterLog([fPrefixLOG fLogStem num2str(scanNo(scan)) fPostfixLOG '.log']); % fPostfix '.log']);
    
    % Read the corresponding Fable g-vector file
    
    peaks_gve_dat = loadGVE([fPrefixGVE fGVEstem num2str(scanNo(scan)) fPostfixGVE '.gve']);
    
    % find g-vectors for grains in the log file, using indices to the g-vector file
    
    clear strn stddevs
    
    igrain = 0;
    
    
    for ig=grainNo(scan):grainNo(scan)
        
        igrain = igrain + 1;
        
        clear gve gveID A B hkl2_
        
        Ngvec = length(filelog(ig).refl(:,1));   % Find the number of g-vectors in the log file
        
        cnt = 0;
        for i=1:Ngvec
            
            if filelog(ig).refl(i,22) <  2.0  % Only include if below some threshold internal angle
                cnt = cnt+1;
                gveID(cnt) = find(peaks_gve_dat(:,9)==filelog(ig).refl(i,3));
                
                % The measured 2-theta, omega and eta for this g-vector
                toe(cnt,:) = filelog(ig).refl(i,[13 16 19]);
                
                hklS(cnt, :) = filelog(ig).refl(i,4:6);
                
            end
        end
        Ngvec = cnt;
        
        % Load in the g-vector and it's norm (or 1/d-spacing)
        gve = peaks_gve_dat(gveID,[1:3 6]);
        
        % The following was a check, to make certain that we did not have
        % Friedel pairs
        % close(1)
        %figure(1)
        %plot(gve(:,2),gve(:,3),'+')
        %hold on
        %plot(-gve(:,2),-gve(:,3),'o');
        
         FriedelPairs = zeros(Ngvec,1);
       
        for i=1:Ngvec
            for j=1:Ngvec
                if i ~= j
                    
                    if isequal(hklS(i,1:3), -hklS(j,1:3))
                        
%                     mag = (dot(a,b)/(norm(a).^2));
%                     if mag>0.998 && mag<1.001

                        FriedelPairs(i) = 1;
                        FriedelPairs(j) = 1;
                    end
                end
            end
        end
   
        % FIX-UP DISTANCE
        % gve = 0.999*gve;
        
        
%         idx = find(FriedelPairs==1);
%         
%         [Ngvec length(idx)]
%         
%         Ngvec = length(idx);
%         gve = gve(idx,:);   
        
        OutSide = zeros(Ngvec,1);
        for i=1:Ngvec
            OutSide(i) = (norm(gve(i,1:3)) > 0.92) & (norm(gve(i,1:3)) < 0.95) ; % 0.795;
        end
        
        idx = find(OutSide ~= 1);
        
        [Ngvec length(idx)];
        
        Ngvec = length(idx);
        gve = gve(idx,:);  
        

        
        
        % Now, the least squares problem
        % A*X = B, or X = A\B        
        
        B = zeros(Ngvec,1);
        % check if the measured g vectors is close to the ideal g vectors
        %figure(2)
        %plot(1:Ngvec,gve(:,4),'+',1:length(g_ideal),g_ideal,'o') % '+": measured g vectors and their norm; 'o': 16 ideal g vectors and their norm based on a and c 
        %hold on
%         for i=1:length(g_ideal)
%             plot([1 Ngvec],[g_ideal(i) g_ideal(i)],'-');
%         end
        %hold off
  
 % pause        

        for i=1:Ngvec
            
            % The delta spacing will be a minimum with corresponding ideal hkl
            gs = norm(gve(i,1:3))-g_ideal;
            [val,midx] = min(abs(gs));
            B(i,1) = -gs(midx)/g_ideal(midx); % B is the measured elastic strain for each reflection (i.e. g vector)
        end
        

        for i=1:Ngvec
            
            lmn = gve(i,1:3)/norm(gve(i,1:3)); % lmn is the g vector in sample coordinate system
            
            dx_term = -( cosd(toe(i,2)) + (sind(toe(i,2))*sind(toe(i,3)))/tand(toe(i,1)/2.0) )/1061117.422;
            dy_term = -( sind(toe(i,2)) + (cosd(toe(i,2))*sind(toe(i,3)))/tand(toe(i,1)/2.0) )/1061117.422;
            
            dx_term = -( cosd(toe(i,2)) + (sind(toe(i,2))*sind(toe(i,3)))/tand(toe(i,1)/1.0) )/999654.8; % new parameter
            dy_term = -( sind(toe(i,2)) + (cosd(toe(i,2))*sind(toe(i,3)))/tand(toe(i,1)/1.0) )/999654.8; % new parameter from Ti_1516.par
            
            %dx_term = -( cosd(toe(i,2)-20) + (sind(toe(i,2)-20)*sind(toe(i,3)))/tand(toe(i,1)/2.0) )/999654.8; % new parameter
            %dy_term = -( sind(toe(i,2)-20) + (cosd(toe(i,2)-20)*sind(toe(i,3)))/tand(toe(i,1)/2.0) )/999654.8; % new parameter from Ti_1516.par
            
            A(i,:) = [lmn(1)^2 lmn(2)^2 lmn(3)^2 2*lmn(1)*lmn(2) 2*lmn(1)*lmn(3) 2*lmn(2)*lmn(3) -dx_term -dy_term];
            
                    
        end
        
        X =A\B; % matlab solve least square problem in one command!
        %figure(4)
        %plot(B-A*X);
        fit_error = B - A*X;
        
        strn(:,igrain) = A\B;
        stddevs(igrain) = std(B-A*X);
        
        % Filter results from the fit, here.
        idx = find(abs(fit_error)<0.004); % .0018); % 0.8e-3); % mark's 0.8e-3);
        %fprintf(1,'%d of %d g-vectors used in the fit for grain %d.\n',length(idx),Ngvec,ig);
        
        B1 = B(idx);
        A1 = A(idx,:);
        
        strn(:,igrain) = A1\B1;
        %figure(5)
        %plot(B1-A1*strn(:,igrain))
        
        stddevs(igrain) = std(B1-A1*strn(:,igrain));
        
        %pause
    end
    
    Rowstrn(scan,:)=strn'; % display strain and position for all grains
    %------------------------------------- below lines are written by Leyun---------------------------------------------------
    strain(1,1)=strn(1); strain(2,2)=strn(2);strain(3,3)=strn(3); strain(1,2)=strn(4);strain(2,1)=strn(4);strain(1,3)=strn(5);strain(3,1)=strn(5);strain(2,3)=strn(6);strain(3,2)=strn(6);
    filelog(ig).U=filelog(ig).U*[0 -1 0;1 0 0;0 0 1]; % 90 deg rotation
    strain_c(:,:)=filelog(ig).U'*strain(:,:,igrain)*filelog(ig).U; % strain tensor in crystal coordinate, filelog(ig).U is the R matrix of the grain
    C=[162.4e3,   92e3,  69e3,  0.0,   0.0,  0.0;
       92e3, 162.2e3,   69e3,  0.0,   0.0,  0.0;
       69e3,  69e3,  181.6e3,  0.0,   0.0,  0.0;
       0.0,     0.0,     0.0,  47.2e3,  0.0,  0.0;
       0.0,     0.0,     0.0,    0.0, 47.2e3, 0.0;
       0.0,     0.0,     0.0,    0.0,   0.0, 35.2e3;]; % stiffness tensor of Ti
    strain_c_vec=[strain_c(1,1), strain_c(2,2),strain_c(3,3),strain_c(2,3)*2,strain_c(3,1)*2,strain_c(1,2)*2]';
    stress_c_vec=C*strain_c_vec; % Stress in the crystal coordinate system using General Hooke's Law
    stress_crystal(scan,:)=stress_c_vec';
    stress_c=[stress_c_vec(1),stress_c_vec(6),stress_c_vec(5);stress_c_vec(6),stress_c_vec(2),stress_c_vec(4);stress_c_vec(5),stress_c_vec(4),stress_c_vec(3)];
    
  
    
    for tn=1:6
        schmid(scan,tn)= T1twin(tn, 1:3)*(filelog(ig).U'*[0 0 0;0 0 0;0 0 1]*filelog(ig).U)*T1twin(tn, 4:6)'; % Schmid factor
        
        twinshear(scan,tn)=mT1twin(1,1,tn)*stress_c(1,1)+mT1twin(1,2,tn)*stress_c(1,2)+mT1twin(1,3,tn)*stress_c(1,3)+mT1twin(2,1,tn)*stress_c(2,1)+mT1twin(2,2,tn)*stress_c(2,2)+mT1twin(2,3,tn)*stress_c(2,3)+mT1twin(3,1,tn)*stress_c(3,1)+mT1twin(3,2,tn)*stress_c(3,2)+mT1twin(3,3,tn)*stress_c(3,3); % resolved shear stress in the twinning system (tn)
        %twinshear(scan,tn)=T1twin(tn, 1:3)*stress_c*T1twin(tn, 4:6)'; % this is actually equivalent to the line above
        twinplane(1:3,tn)=filelog(ig).U*T1twin(tn, 1:3)'; % twin plane normal in the sample coordinate
        twindr(1:3,tn)=filelog(ig).U*T1twin(tn, 4:6)'; % twin direction in the sample coordinate
    end
    if tn==6
        
        twinplane;
        twindr; % output twin plane and twin direction for all six twinning variants
    end
    stress=filelog(ig).U*stress_c*filelog(ig).U'; % stress in the sample coordinate system
    stress_vec(scan,:)=[stress(1,1),stress(2,2),stress(3,3),stress(1,2),stress(1,3),stress(2,3)]'; % stress in the sample coordinate system in vector form
        
    % --------------------------------------------------end of Leyun's coding--------------------------------------------
    EulerAng(scan,1)=log_data(scan).euler(1);
    EulerAng(scan,2)=log_data(scan).euler(2);
    EulerAng(scan,3)=log_data(scan).euler(3)+30;
    
    GrainArray(scan,1)=scan;
    GrainArray(scan,2:3)=Rowstrn(scan,7:8);
    GrainArray(scan,4:9)=stress_vec(scan,:);
    GrainArray(scan,10:12)=EulerAng(scan,:);


end









end