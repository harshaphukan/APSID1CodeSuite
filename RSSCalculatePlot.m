
clc
clear
close all
load('Layer10ParentGrainNum.mat')
load('L10GrainInfo.mat')
load('ID1StrainInfo.mat')
load('S.mat')
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
Strain=ID1StrainInfo(:,2);

I1=find(ID1StrainInfo(:,1)~=282 & ID1StrainInfo(:,1)~=326);
Strain=Strain(I1);
RSS=RSS(I1,:);

%%
%Plotting the RSS for the different twin variants
Orange=[255,69,0]/256; % Orange Color RGB
Turq=[64,224,208]/256; % Turquiose Color RGB
%sz=140; % Not necessary
figure
plot(Strain,RSS(:,1),'bd')
hold on
plot(Strain,RSS(:,2),'ms')

hold on
plot(Strain,RSS(:,3),'^','MarkerEdgeColor',Turq)

hold on
plot(Strain,RSS(:,4),'x','MarkerEdgeColor',Orange)
hold on
%plot(ID1StrainInfo(:,2),RSS(:,5),'g+')
scatter(Strain,RSS(:,5),'ro','Filled')
hold on
plot(Strain,RSS(:,6),'ko')

hold off
xlabel('Bulk Strain/%')
ylabel('Resolved Shear Stress/MPa')
ylim([-100,400])
%title('Layer11:Resolved shear stress Evolution')
grid on
legend('Variant 1','Variant 2','Variant 3','Variant 4','Variant 5','Variant 6',...
      'Location','northwest')

print('-f1','Layer10RSSPlotNoTitle','-dpdf')

%GrainInfo(cellfun('isempty', GrainInfo)) = [];
