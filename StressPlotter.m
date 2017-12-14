clc
%clear
close all
%D=xlsread('Layer01SigmaZZGrain12.xlsx');

load('L2ParentL3GrainEuclid.mat')
load('Layer3FilesGrainInfo.mat')
load('ID1StrainInfo.mat')
load('S.mat')
%%
L=length(GrainNum);
STens=zeros(L,6);
for i=1:L
    if GrainNum(i)~=0
    STens(i,:)=GrainInfo{i}(GrainNum(i),4:9);
    else
        STens(i,:)=zeros(1,6);
    end
end
%%
Strain=ID1StrainInfo(:,2);

%%Plotting the Stress
%{
Strain(14)=[];
STens(14,:)=[];
%}
figure
scatter(Strain,STens(:,3),'bo','filled')

xlabel('Bulk Strain/%')
ylabel('{\sigma_{zz}}/MPa')
ylim([-500,600])
%title('Layer11:Resolved shear stress Evolution')
grid on

print('-f1','Layer2SigmaZZPlotNoTitle','-dpdf')

%GrainInfo(cellfun('isempty', GrainInfo)) = [];
