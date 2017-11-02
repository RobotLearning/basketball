%% Analyze real robot joint data

% load joints data
clc; clear; close all;

filename = '~/basketball/data/joints.txt';
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

%% Draw the cartesian desired and actual positions
cart_filename = '~/basketball/data/cartesian.txt';
C = dlmread(cart_filename);
figure('Name','Cartesian positions');
downsample = 5;
vec = 1:5:size(C,1);
left_des_pos = C(vec,1:3);
right_des_pos = C(vec,4:6);
left_act_pos = C(vec,13:15);
right_act_pos = C(vec,16:18);
hold on;
grid on;
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
scatter3(left_des_pos(:,1),left_des_pos(:,2),left_des_pos(:,3),'r');
scatter3(left_act_pos(:,1),left_act_pos(:,2),left_act_pos(:,3),'k');
scatter3(right_des_pos(:,1),right_des_pos(:,2),right_des_pos(:,3),'r');
scatter3(right_act_pos(:,1),right_act_pos(:,2),right_act_pos(:,3),'b');

% compare with SL computed kinematics
% cart_sl_filename = '~/basketball/data/cartesian_SL.txt';
% CSL = dlmread(cart_sl_filename);
% left_act_sl = CSL(vec,1:3);
% right_act_sl = CSL(vec,7:9);
% scatter3(left_act_sl(:,1),left_act_sl(:,2),left_act_sl(:,3),'g');
% scatter3(right_act_sl(:,1),right_act_sl(:,2),right_act_sl(:,3),'g');

%% Draw the ball and initial robot loc
string_len = 1.0;
ball_radius = 0.1213;
ball_loc = [0.0, 0.5, 1.5-string_len-ball_radius];
basketball_color = [207,83,0]/256;
numPoints = 100;
[ballMeshX,ballMeshY,ballMeshZ] = sphere(numPoints);
ballSurfX = ball_loc(1) + ball_radius * ballMeshX;
ballSurfY = ball_loc(2) + ball_radius * ballMeshY;
ballSurfZ = ball_loc(3) + ball_radius * ballMeshZ;
h = surf(ballSurfX,ballSurfY,ballSurfZ);
set(h,'FaceColor',basketball_color,'FaceAlpha',0.8,'EdgeAlpha',0.05);

robot_left_init = left_act_pos(1,:);
robot_right_init = right_act_pos(1,:);
numPoints = 10;
[robotInitMeshX,robotInitMeshY,robotInitMeshZ] = sphere(numPoints);
robotInitSurfX = robot_left_init(1) + 0.02 * robotInitMeshX;
robotInitSurfY = robot_left_init(2) + 0.02 * robotInitMeshY;
robotInitSurfZ = robot_left_init(3) + 0.02 * robotInitMeshZ;
h = surf(robotInitSurfX,robotInitSurfY,robotInitSurfZ);
robotInitSurfX = robot_right_init(1) + 0.02 * robotInitMeshX;
robotInitSurfY = robot_right_init(2) + 0.02 * robotInitMeshY;
robotInitSurfZ = robot_right_init(3) + 0.02 * robotInitMeshZ;
h = surf(robotInitSurfX,robotInitSurfY,robotInitSurfZ);

legend('left des', 'left act', 'right des', 'right act', ...
       'left init pos', 'right init pos');