%% Estimate calibration matrix

% load the data
clc; clear; close all;

train_file = '~/basketball/data/blobs_tabletennis.txt';
test_file = '~/basketball/data/blobs_tabletennis_2.txt';
yellow_ball_file = '~/basketball/data/blobs_yellowball.txt';
M = dlmread(train_file);
N = size(M,1);
% get the robot left arm cartesian endeffector positions
R = M(:,2:4)';
Rbar = [R; ones(1,N)];
% get camera basketball center of mass estimates
C = M(:,5:end)';
Cbar = [C; ones(1,N)];

% estimate with constrained nonlinear least squares
%{
M1 = 1000*ones(2,3);
M2 = 1000*ones(2,3);
o1 = -M1 * [0.05; 0.00; 0.5];% from camera 1 to robot base (origin)
o2 = -M2 * [-0.05; 0.00; 0.5]; % from camera 2 to robot base
theta = zeros(3,1); % Euler angles
d = [0.2; 0.3; 0.0];
%x0 = rand(24,1);
x0 = [M1(:); M2(:); o1; o2; theta; d];
fun = @(x) fit_calibration_matrix(Cbar,Rbar,x);
options = LMFnlsq('default');
options = LMFnlsq(options,'Display',10,'MaxIter',1000);
x = LMFnlsq(fun,x0,options);
[M,D] = form_calibration_matrices(x);
Mat = M*D;
%}
% estimate with least squares
%Mat = C / R;
Mat = R / Cbar;

%% Check the estimation results for training
Rest = Mat * Cbar;
%Rest = Mat \ Cbar;
%Rest = Mat \ C;
dt = 0.002;
t = dt * (1:size(Rest,2));
rms_train = sqrt((norm(Rest(1:3,:) - R, 'fro')^2) / (size(Rest,2)))

figure('Name','Training data');
subplot(3,1,1);
plot(t, Rest(1,:),'-b',t, R(1,:),'--r');
ylabel('Robot hand pose X');
subplot(3,1,2);
plot(t, Rest(2,:),'-b',t, R(2,:),'--r');
ylabel('Robot hand pose Y');
subplot(3,1,3);
plot(t, Rest(3,:),'-b',t, R(3,:),'--r');
ylabel('Robot hand pose Z');

%{
figure('Name','Training data');
subplot(2,2,1);
plot(t, Cest(1,:),'-b',t, C(1,:),'--r');
ylabel('Camera 1-X');
subplot(2,2,2);
plot(t, Cest(2,:),'-b',t, C(2,:),'--r');
ylabel('Camera 1-Y');
subplot(2,2,3);
plot(t, Cest(3,:),'-b',t, C(3,:),'--r');
ylabel('Camera 2-X');
subplot(2,2,4);
plot(t, Cest(4,:),'-b',t, C(4,:),'--r');
ylabel('Camera 2-Y');
%}

%% Check results for test data
Mtest = dlmread(test_file);
Ntest = size(Mtest,1);
Rtest = Mtest(:,2:4)';
Rtest_bar = [Rtest; ones(1,Ntest)];
Ctest = Mtest(:,5:end)';
Ctest_bar = [Ctest ; ones(1,Ntest)];
t_test = dt * (1:size(Ctest,2));
Rest_test = Mat * Ctest_bar;
%Rest_test = Mat \ Ctest_bar;
%Rest_test = Mat \ Ctest;

rms_test = sqrt((norm(Rest_test(1:3,:) - Rtest, 'fro')^2) / size(Rtest,2))

Mtest2 = dlmread(yellow_ball_file);
Ntest2 = size(Mtest2,1);
Cyellow_test = Mtest2(:,5:end)';
Cyellow_test_bar = [Cyellow_test; ones(1,Ntest2)];
Ytest = Mat * Cyellow_test_bar;
t_test_yellow = dt * (1:Ntest2);
figure('Name','Test data - yellow ball');
subplot(3,1,1);
plot(t_test_yellow, Ytest(1,:),'-b');
ylabel('Yellow ball X');
subplot(3,1,2);
plot(t_test_yellow, Ytest(2,:),'-b');
ylabel('Yellow ball Y');
subplot(3,1,3);
plot(t_test_yellow, Ytest(3,:),'-b');
ylabel('Yellow ball Z');

figure('Name','Test data - yellow ball cartesian');
plot3(Ytest(1,:),Ytest(2,:),Ytest(3,:));
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;
grid on;