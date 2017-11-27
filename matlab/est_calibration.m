%% Estimate calibration matrix

% load the data
clc; clear; close all;

train_file = '~/basketball/data/blobs1.txt';
test_file = '~/basketball/data/blobs2.txt';
M = dlmread(train_file);
N = size(M,1);
% get the robot left arm cartesian endeffector positions
R = M(:,2:4)';
Rbar = [R; ones(1,N)];
% get camera basketball center of mass estimates
C = M(:,5:end)';
Cbar = [C; ones(1,N)];

% estimate with constrained nonlinear least squares
m0 = rand(3,5);
x0 = m0(:);
fun = @(x) fit_calibration_matrix(R,Cbar,x);
options = LMFnlsq('default');
options = LMFnlsq(options,'Display',1,'MaxIter',100);
x = LMFnlsq(fun,x0,options);
Mat = reshape(x,3,5);

%% Check the estimation results for training
Rest = Mat * Cbar;
dt = 0.002;
t = dt * (1:size(Rest,2));
% 
% figure('Name','Training data');
% subplot(2,2,1);
% plot(t, Cest(1,:),'-b',t, C(1,:),'--r');
% ylabel('Camera 1-X');
% subplot(2,2,2);
% plot(t, Cest(2,:),'-b',t, C(2,:),'--r');
% ylabel('Camera 1-Y');
% subplot(2,2,3);
% plot(t, Cest(3,:),'-b',t, C(3,:),'--r');
% ylabel('Camera 2-X');
% subplot(2,2,4);
% plot(t, Cest(4,:),'-b',t, C(4,:),'--r');
% ylabel('Camera 2-Y');

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

%% Check results for test data
Mtest = dlmread(test_file);
Ntest = size(Mtest,1);
Rtest = Mtest(:,2:4)';
Rtest = [Rtest; ones(1,Ntest)];
Ctest = Mtest(:,5:end)';
Ctest_bar = [Ctest ; ones(1,Ntest)];
t_test = dt * (1:size(Ctest,2));
Rest_test = Mat * Ctest_bar;


figure('Name','Test data');
subplot(3,1,1);
plot(t_test, Rtest(1,:),'-b',t_test, Rest_test(1,:),'--r');
ylabel('Robot hand pose X');
subplot(3,1,2);
plot(t_test, Rtest(2,:),'-b',t_test, Rest_test(2,:),'--r');
ylabel('Robot hand pose Y');
subplot(3,1,3);
plot(t_test, Rtest(3,:),'-b',t_test, Rest_test(3,:),'--r');
ylabel('Robot hand pose Z');