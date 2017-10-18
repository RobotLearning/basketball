%% analyze real robot data

% load joints data
clc; clear; close all;

filename = '~/basketball/joints.txt';
M = dlmread(filename);
t = 1:size(M,1);
figure('Name','LEFT ARM');
q_des = M(:,1:7);
qd_des = M(:,15:21);
q_act = M(:,29:35);
qd_act = M(:,43:49);
for i = 1:7
    subplot(7,2,(i-1)*2 + 1);
    plot(t,q_des(:,i),'r',t,q_act(:,i),'b');
    subplot(7,2,(i-1)*2 + 2);
    plot(t,qd_des(:,i),'r',t,qd_act(:,i),'b');    
end
figure('Name','RIGHT ARM');
q_des = M(:,8:14);
qd_des = M(:,22:28);
q_act = M(:,36:42);
qd_act = M(:,50:56);
for i = 1:7
    subplot(7,2,(i-1)*2 + 1);
    plot(t,q_des(:,i),'r',t,q_act(:,i),'b');
    subplot(7,2,(i-1)*2 + 2);
    plot(t,qd_des(:,i),'r',t,qd_act(:,i),'b');    
end