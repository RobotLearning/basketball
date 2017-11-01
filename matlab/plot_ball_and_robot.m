%% Plotting ball and robot

clc; clear; close all;

idx = 1:650;
P = load('balls_pos.txt');
%V = load('balls_vel.txt');
R = load('robot_cart.txt');

basketball_color = [207,83,0]/256;
ball_radius = 0.1213;
numPoints = 100;
[ballMeshX,ballMeshY,ballMeshZ] = sphere(numPoints);
ballSurfX = P(1,idx(end)) + ball_radius * ballMeshX;
ballSurfY = P(2,idx(end)) + ball_radius * ballMeshY;
ballSurfZ = P(3,idx(end)) + ball_radius * ballMeshZ;

figure;
scatter3(P(1,idx),P(2,idx),P(3,idx),'b');
hold on;
axis equal; 
%scatter3(M(1,v(end)),M(2,v(end)),M(3,v(end)),1000,basketball_color,'fill');
h = surf(ballSurfX,ballSurfY,ballSurfZ);
set(h,'FaceColor',basketball_color,'FaceAlpha',1,'EdgeAlpha',0);
hold on;
scatter3(R(1,idx),R(2,idx),R(3,idx),'r')
hold on;
scatter3(R(4,idx),R(5,idx),R(6,idx),'k');
legend('balls past','last ball pos','left arm', 'right arm');

% figure;
% scatter3(V(1,idx),V(2,idx),V(3,idx),'b');
% xlabel('x-vel');
% ylabel('y-vel');
% zlabel('z-vel');
% hold on;
% dt = 0.002;
% vel_diff = diff(P')'/dt;
% scatter3(vel_diff(1,idx),vel_diff(2,idx),vel_diff(3,idx),'r');