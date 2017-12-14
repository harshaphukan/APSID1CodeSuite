%% Script to Plot Twinning Events as a Function of the Reported Resolved Shear Stress Associated with the Activated Twin Variant.

clc
clear
close all

A=importdata('RSSMPrimeDataID1Experiment_July2013.xlsx'); % Load .xlsx file with Twin/SS Data
SSys=A.textdata(:,end);
SSys=SSys(2:end);
UniqueSlipSys=unique(SSys); % Unique Types of Slip Systems Observed

%SSys=C(2:end,end);
%UniqueSlipSyS=unique(SSys); % Unique Types of Slip Systems Observed

LayerNum=A.data(:,1); % Layer Number
%%
SFNeighbor=A.data(:,2); % Schmid Factor of Neighbor Grain with High m-prime correlation
SFRank=A.data(:,3); % Schmid Factor Rank of Activated T1 Twin Variant [6 Highest; 1 Lowest]
RSS=A.data(:,4); % Resolved Shear Stress/MPa
mprime=A.data(:,5);

PlausVec=A.data(:,6); % Plausibility of Slip/Twin Transfer Across Grain Bdy [Plausible=1; Not Plausible=0]

RVec=1:6; % Ranking of Twin Variants by Schmid Factor: 1 is Lowest: 6 is Highest
RVec=50*RVec; % Multiply Schmid Vector Rank vector by a weighted value!
XRSS=linspace(0,1,200);
MaxRSS=max(RSS)*ones(length(XRSS));
MinRSS=min(RSS)*ones(length(XRSS));

% Color Convention on the basis of Slip System: Blue for basal, Red for
% Prism, Green for pyr-a and Gold for pyr <c+a>
F=0; % Initialize The filled/open circle string containing variable 

 ColorVec=zeros(1,3,5);
 ColorVec(:,:,1)=[0 0 139]/255;
ColorVec(:,:,2)=[255,215,0]/255;
ColorVec(:,:,3)=[128 0 0]/255;
ColorVec(:,:,4)=[255 0 0]/255;
ColorVec(:,:,5)=[124 252 0]/255;
 

C=zeros(1,3); % Variable that holds color string
figure(1)


for ii=1:length(RSS)
     
    % Assign Colors Based on Type of Slip System Observed
     if strcmp(SSys(ii),'basal')==1
         C=[0 0 139]/255;% Blue for basal
        
     elseif strcmp(SSys(ii),'pyr<c+a>')==1 || strcmp(SSys(ii),'pyr <c+a>')==1 
        C=[255,215,0]/255; % Gold for pyr <c+a>
    elseif strcmp(SSys(ii),'T1')==1
        C=[128 0 0]/255; % Dark Red for T1 System
     elseif strcmp(SSys(ii),'prism')==1
        C=[255 0 0]/255; % Red for Prism Slip System
    else
        C=[124 252 0]/255; % Lawn Green for Pyramidal<a>
     end
         if PlausVec(ii)==1
             F='filled';
        p1=scatter(mprime(ii),RSS(ii),SFRank(ii)*50,C,F);
        
         else
             
         p1=scatter(mprime(ii),RSS(ii),SFRank(ii)*50,C);
         
        
            
        end
        hold on
        grid on
    
  
end
 k1=plot(XRSS,MaxRSS+10,'b','LineWidth',2.0); 
 k2=plot(XRSS,MinRSS-10,'r','LineWidth',2.0);
h = zeros(5, 1);
h(1) = scatter(NaN,NaN,50,ColorVec(:,:,1),'filled');
h(2) = scatter(NaN,NaN,50,ColorVec(:,:,2),'filled');
h(3) = scatter(NaN,NaN,50,ColorVec(:,:,3),'filled');
h(4) = scatter(NaN,NaN,50,ColorVec(:,:,4),'filled');
h(5) = scatter(NaN,NaN,50,ColorVec(:,:,5),'filled');

legend(h, 'Basal','Pyramidal <c+a>','T1','Prism','Pyramidal <a>','Location','southwest'); 

       

P={'Slip System Interaction of Parent Grain with Neighbors'};
axis([0.4 1 1 350])
whitebg('w')
xlabel('Slip Transfer Parameter(m'')','FontSize',10)
ylabel('Resolved Shear Stress/MPa','FontSize',10)
title(P,'FontSize',12)
xt=get(gca,'XTick');
%yt=get(gca,'YTick');
set(gca,'FontSize',12,'FontWeight','bold');
set(gca,'ticklength',2.0*get(gca,'ticklength'))
set(gca,'YTick',0:100:300);
box on


print('-f1','AllTwinRSSMPrimePlot','-dpdf')




