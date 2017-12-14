clc
FilePath='/Users/harshaphukan/Documents/LaTeXProjects/PaperManuscriptDraft2017/';
filename='3DRecons_052114.seeds';
fid=fopen(fullfile(FilePath,filename),'r');
if(fid == -1)
    beep;
    error('Cannot open file:\n  %s\n', fileName);
end


textscan(fid, '%*[^\n]', 5); % Skip headers


  [M,Count] = fscanf(fid,'%f');

  OutArray=zeros(Count/6,6);
  i=1;
  while i<=Count
      for j=1:Count/6
  OutArray(j,:)=M(i:i+5);
 i=i+6;
 OutArray(j,6)=OutArray(j,6)+30;
      end
  end
  
 
fclose(fid);
%%
Euler9=L9_98A(78,10:12);

E9=[84.7696   15.2082   95.5605];
MisOr=zeros(length(OutArray),1);
for i=1:length(OutArray)
    MisOr(i)=HexMisOr(E9,OutArray(i,4:6));
end

IV=find(MisOr<5);

 % Candidates that meet MisOrientation Criteria
    
for j=1:length(L9_98A)
    MisOr2(j)=HexMisOr(E9,L9_98A(j,10:12));
end
IV1=find(MisOr2<5)
