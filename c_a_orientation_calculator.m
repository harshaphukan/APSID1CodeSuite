%% This script calculates the c and a orientations between two grains: Requires the data structure from .log file (log_data) as input
function [C_ang,A1_ang,A2_ang,A3_ang,A4_ang]=c_a_orientation_calculator(log_data)
GN1=input('Enter number of first grain: ','s');
GN2=input('Enter number of second grain: ','s');
GA=str2num(GN1);
GB=str2num(GN2);
euler_A=log_data(GA).euler;
euler_B=log_data(GB).euler;
%Construction of Rotation matrices
R_phi1_A=[cosd(euler_A(1)),-sind(euler_A(1)),0;sind(euler_A(1)),cosd(euler_A(1)),0;0,0,1];
R_phi_A=[1,0,0;0,cosd(euler_A(2)),-sind(euler_A(2));0,sind(euler_A(2)),cosd(euler_A(2))];
R_phi2_A=[cosd(euler_A(3)+30),-sind(euler_A(3)+30),0;sind(euler_A(3)+30),cosd(euler_A(3)+30),0;0,0,1];%With 30-degree correction for phi2
R_A=R_phi1_A*R_phi_A*R_phi2_A; % Rotation matrix for grain A
R_phi1_B=[cosd(euler_B(1)),-sind(euler_B(1)),0;sind(euler_B(1)),cosd(euler_B(1)),0;0,0,1];
R_phi_B=[1,0,0;0,cosd(euler_B(2)),-sind(euler_B(2));0,sind(euler_B(2)),cosd(euler_B(2))];
R_phi2_B=[cosd(euler_B(3)+30),-sind(euler_B(3)+30),0;sind(euler_B(3)+30),cosd(euler_B(3)+30),0;0,0,1];%With 30-degree correction for phi2
R_B=R_phi1_B*R_phi_B*R_phi2_B; % Rotation matrix for grain B
%Calculation of orientation vectors
c_A=R_A*[0 0 1]';
c_B=R_B*[0 0 1]';

a1_A=R_A*[1;0;0];
a2_A=R_A*[cosd(120);sind(120);0];

a1_B=R_B*[1;0;0];
a2_B=R_B*[cosd(120);sind(120);0];

C_ang=acosd(dot(c_A,c_B));

A1_ang=acosd(dot(a1_A,a1_B));
A2_ang=acosd(dot(a1_A,a2_B));
A3_ang=acosd(dot(a2_A,a1_B));
A4_ang=acosd(dot(a2_A,a2_B));

end

