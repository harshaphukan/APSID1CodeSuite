% Script to generate Grain numbers of the same grain in successive load steps..checks for the smallest misorientation 
% with the compared grain
clc
clear
FilePath=input('Enter File Path\n','s');
LF=load('LogFileListLayer3.mat');
LogFile=LF.LogFile;
%save('LogFileListL3P2.mat','LogFile');
%%      
for i=1:length(LogFile)
GrainInfo{i}=LoadDataSameDir(FilePath,LogFile{i});

end
save('Layer3FilesGrainInfo.mat','GrainInfo');
%%
load('Layer3FilesGrainInfo.mat');
%%
Euler1=[-4.0732,331.1366,80.4394];
D=[-129.8129  360.4559];
%P1=GrainInfo{20}(12,2:3);
%%
GrainNum=zeros(1,length(LogFile));
for i=1:length(LogFile)
for j=1:length(GrainInfo{i})
MisOr{i}(j)=HexMisOr(Euler1,GrainInfo{i}(j,10:12));
%GrainNum(i)=find(MisOr{i}==min(MisOr{i}));
Dist{i}(j)=Euclidean(D,GrainInfo{i}(j,2:3));
IOr=MisOr{i}<=5; IDist=Dist{i}<130;
I=IOr.*IDist;
F1=find(I==1);
if isempty(F1)
GrainNum(i)=0;%find(MisOr{i}<=5 && Dist{i}<=140);
else
    GrainNum(i)=F1(1);
end
end
end
%%
save('L2ParentL3GrainEuclid.mat','GrainNum') 

%%
%{
ID1StrainInfo(:,1)=[90,101,112,127,138,149,160,171,183,194,205,227,238,249,260,...
                    271,282,293,304,315,326,337,348,359,370,381,392,403,416,427,438,...
                    449,460,471,482,493,504];

ID1StrainInfo(:,2)=[0,0.03,0.03,0.07,0.14,0.22,0.23,0.25,0.31,0.36,0.37,0.46,0.53,0.64,...
                  0.72,0.85,0.98,1.11,1.28,1.48,1.56,1.76,1.89,2.04,2.17,2.42,2.50,2.60,...
                  2.73,2.93,2.91,2.88,2.84,2.8,2.73,2.66,2.59];  
save('ID1StrainInfo.mat','ID1StrainInfo')
%}