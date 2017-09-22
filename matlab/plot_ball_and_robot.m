%% Plotting ball and robot

clc; clear; close all;

v = 1:1000;
M = load('balls_pred.txt');
R = load('robot_cart.txt');
orange = [255,153,0]/256;

figure;
scatter3(M(1,v),M(2,v),M(3,v),'b');
hold on;
scatter3(M(1,v(end)),M(2,v(end)),M(3,v(end)),100,orange,'fill');
hold on;
scatter3(R(1,v),R(2,v),R(3,v),'r')
hold on;
scatter3(R(4,v),R(5,v),R(6,v),'k');
legend('balls past','last ball pos','right arm', 'left arm');