function [x,y] = make_training_data()

x = [];
y = [];

% cluster 1
mu = [2,3];
sigma = [1,1.5;1.5,3];
r = mvnrnd(mu,sigma,100);


plot(r(:,1),r(:,2), 'bo', 'MarkerFaceColor', 'b');
x = [x; r(:,1)]; y = [y; r(:,2)];

% cluster 2
mu = [8,6];
sigma = [2 ,  -1.5
        -1.5,  3];
r = mvnrnd(mu,sigma,100);
plot(r(:,1),r(:,2), 'ro', 'MarkerFaceColor', 'r');
x = [x; r(:,1)]; y = [y; r(:,2)];

% cluster 3
mu = [13,4];
sigma = [1.05 ,  0.75
        0.75,  1.05];
r = mvnrnd(mu,sigma,100);
plot(r(:,1),r(:,2), 'go', 'MarkerFaceColor', 'g');
x = [x; r(:,1)]; y = [y; r(:,2)];

% cluster 4
mu = [20,-3];
sigma = [1,-1.5;-1.5,3];
r = mvnrnd(mu,sigma,100);
plot(r(:,1),r(:,2), 'bo', 'MarkerFaceColor', 'b');
x = [x; r(:,1)]; y = [y; r(:,2)];

% cluster 5
mu = [45, -20];
sigma = [1.05 ,  0.75
        0.75,  1.05];
r = mvnrnd(mu,sigma,100);
plot(r(:,1),r(:,2), 'go', 'MarkerFaceColor', 'g');
x = [x; r(:,1)]; y = [y; r(:,2)];


% cluster 6
mu = [30, 10];
sigma = [1.05 ,  0.75
        0.75,  1.05];
r = mvnrnd(mu,sigma,100);
plot(r(:,1),r(:,2), 'go', 'MarkerFaceColor', 'g');
x = [x; r(:,1)]; y = [y; r(:,2)];




%x = [r1(:,1);r2(:,1);r3(:,1)];
%y = [r1(:,2);r2(:,2);r3(:,2)];

%plot(x,y, 'k+');

