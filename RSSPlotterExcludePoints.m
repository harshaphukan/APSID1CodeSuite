clc
clear
load('Layer10ParentGrainNum.mat')
load('L10GrainInfo.mat')
load('ID1StrainInfo.mat')
load('S.mat')
%%
%{
GrainNum(17)=[];%GrainNum(20)=[];
ID1StrainInfo(17,:)=[]; %ID1StrainInfo(20,:)=[];
GrainInfo{17}=[]; %GrainInfo{17}=[];
GrainInfo(cellfun('isempty', GrainInfo)) = [];
%}
%%
L=length(GrainNum);
STens=zeros(L,6);
for i=1:L
    STens(i,:)=GrainInfo{i}(GrainNum(i),4:9);
end


%%
mT1=zeros(6,3);
nT1=zeros(6,3);

for i=1:6
mT1(i,:)=S(i+33).f2; % Shear Direction unit vectors for T1 variants 34-39
nT1(i,:)=S(i+33).f3;% Plane Normal unit vectors for T1 variants 34-39
end


RSS=zeros(L,6); % RSS of 6 T1 Variants 34--->39
for i=1:L
    for j=1:6
RSS(i,j)=RSSEvaluator(mT1(j,:),nT1(j,:),STens(i,:));
    end
end

%%
% Remove Selected outlier points
RSS(17,:)=[]; RSS(20,:)=[];
Strain=ID1StrainInfo(:,2);
Strain(17)=[]; Strain(20)=[];


%%
%Plotting the RSS for the different twin variants

figure
plot(Strain,RSS(:,1),'b^')
hold on
plot(Strain,RSS(:,2),'mh')
hold on
plot(Strain,RSS(:,3),'gs')
hold on
plot(Strain,RSS(:,4),'kv')
hold on
scatter(Strain,RSS(:,5),'ro','filled')
hold on
plot(Strain,RSS(:,6),'cd')
hold off
xlabel('Bulk Strain/%')
ylabel('Resolved Shear Stress/MPa')
ylim([-100,400])
%title('Layer11:Resolved shear stress Evolution')
grid on
legend('Variant 1','Variant 2','Variant 3','Variant 4','Variant 5','Variant 6',...
      'Location','northwest')

print('-f1','Layer10RSSPlotNoTitle','-dpdf')


