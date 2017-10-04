%% Plotting ball and robot

clc; clear; close all;

v = 1:600;
M = load('balls_pred.txt');
R = load('robot_cart.txt');

basketball_color = [207,83,0]/256;
ball_radius = 0.1213;
numPoints = 100;
[ballMeshX,ballMeshY,ballMeshZ] = sphere(numPoints);
ballSurfX = M(1,v(end)) + ball_radius * ballMeshX;
ballSurfY = M(2,v(end)) + ball_radius * ballMeshY;
ballSurfZ = M(3,v(end)) + ball_radius * ballMeshZ;

figure;
scatter3(M(1,v),M(2,v),M(3,v),'b');
hold on;
axis equal; 
%scatter3(M(1,v(end)),M(2,v(end)),M(3,v(end)),1000,basketball_color,'fill');
h = surf(ballSurfX,ballSurfY,ballSurfZ);
set(h,'FaceColor',basketball_color,'FaceAlpha',1,'EdgeAlpha',0);
hold on;
scatter3(R(1,v),R(2,v),R(3,v),'r')
hold on;
scatter3(R(4,v),R(5,v),R(6,v),'k');
legend('balls past','last ball pos','left arm', 'right arm');