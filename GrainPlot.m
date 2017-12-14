function GrainPlot(Array,LayerNum)
%% This part rotates the COM Map by 25 degrees (CCW)
[R,~]=size(Array);
X=zeros(R,1);
Y=zeros(R,1);
XY25=zeros(R,2);
theta=deg2rad(25); % Rotate COM Positions by 30 degrees CCW in the XY Plane
%Define Orientation Matrix in 2D
RotMat=[cos(theta),-sin(theta);sin(theta),cos(theta)];

for i=1:R
  X(i)=Array(i,2);
  Y(i)=Array(i,3);
  AxRot=RotMat*[X(i),Y(i)]';
  XY25(i,1)=AxRot(1);
  XY25(i,2)=AxRot(2);

end



%% The portion is to plot the spatial COMs of the grains
figure(1)

plot(XY25(:,1),XY25(:,2),'b.','Markersize',15)
%plot(Rowstrn(:,7),Rowstrn(:,8),'b.','Markersize',15)
for i=1:R
   text(XY25(i,1)+3, XY25(i,2)+10,num2str(Array(i,1)),'Fontsize',8);
end

offset=50.0;
xlim([min(XY25(:,1))-offset,max(XY25(:,1))+offset]);
ylim([min(XY25(:,2))-offset,max(XY25(:,2))+offset]);

%axis([-200,0,400,900])
axis square
xlabel('X (micron)','Fontsize',18)
ylabel('Y (micron)','Fontsize',18)
title(['Grain Map' LayerNum],'Fontsize',25)
end
