%% Function to Plot Relative Positions of COMS for two grains in 2D
function Distance=GrainPlotter(COM1,COM2,GN1,GN2)
Distance=sqrt((COM1(1)-COM2(1))^2+(COM1(2)-COM2(2))^2);
figure
plot(COM1(1),COM1(2),'b.','Markersize',15)
text(COM1(1)+3, COM1(2)+10,num2str(GN1),'Fontsize',8);
hold on
plot(COM2(1),COM2(2),'b.','Markersize',15)
text(COM2(1)+3, COM2(2)+10,num2str(GN2),'Fontsize',8);
hold off
title('Relative Positions of Two Grains in 2D Space')
xlabel('X')
ylabel('Y')
xlim([-600,600])
ylim([-600,600])
grid 'on'

end