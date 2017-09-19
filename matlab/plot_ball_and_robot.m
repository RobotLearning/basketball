%% Plotting ball and robot

clc; clear; close all;

v = 1:500;
t = 0.002:0.002:2;
M = load('balls_pred.txt');
R = load('robot_cart.txt');

figure;
% scatter3(M(1,v),M(2,v),M(3,v));
% hold on;
scatter3(R(1,v),R(2,v),R(3,v),'r')