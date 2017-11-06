%% Testing LQR for catching a ball

clc; clear; close all; 

n = 2; % dim_x
m = 2; % dim_u
T = 1.0; % final time
N = 20;
dt = T/N; % discretization
t = dt*(0:N);
Ax = [zeros(n), eye(n); randn(n,2*n)]; % for the hand
Bx = [zeros(n,m); eye(n,m)];
Ab = [zeros(n), eye(n); randn(n,2*n)]; % for the ball
lambda = 1e-3;   
R = lambda * eye(m); 
Q = zeros(2*2*n);
Tr = [eye(n); zeros(n); -eye(n); zeros(n)];
q = 100;
Qf = q * Tr*Tr';

A = [Ax, zeros(2*n); zeros(2*n), Ab];
B = [Bx;zeros(2*n,m)];
C = eye(2*2*n);
b0 = [ones(n,1); zeros(n,1)];
x0 = zeros(2*n,1);

lqr = LQR(Q,R,Qf,A,B,C,N,dt,1);
K = lqr.computeFinHorizonLTI();

x = zeros(2*n,N+1);
b = zeros(2*n,N+1);
x(:,1) = x0;
b(:,1) = b0;
for i = 1:N
    b(:,i+1) = Ab*b(:,i);
    x(:,i+1) = Ax*x(:,i) + Bx*K(:,:,i)*[x(:,i);b(:,i)];
end

plot(t,x(1:n,:),'-',t,b(1:n,:),'--');
legend('robotx','roboty','ballx','bally');
    