%% Function that is used to fit matrices in calibration

% x are the parameters
% y is the camera output
% residual is the vector of differences between y and predicted output
function res = fit_calibration_matrix(y,data,x)

M1 = reshape(x(1:9),3,3);
M2 = reshape(x(10:18),3,3);
o1 = x(19:21);
o2 = x(22:24);
theta = x(25:27); % Euler angles
d = x(28:30); % from camera to origin
M = [M1, -M1*o1; M2, -M2*o2; zeros(1,3), 1];
c1 = cos(theta(1));
s1 = sin(theta(1));
c2 = cos(theta(2));
s2 = sin(theta(2));
c3 = cos(theta(3));
s3 = sin(theta(3));
R = [c1*c2*c3 - s1*s3, -c1*c2*s3 - s1*c3, c1*s2; ...
     s1*c2*c3 + c1*s3, -s1*c2*s3 + c1*c3, s1*s2; ...
     -s2*c3, s2*s3, c2];
D = [R, d; zeros(1,3), 1];
res = y - M*D*data;
% m1 = x(1:n*n);
% m2 = x(n*n+1:end);
% M1 = reshape(m1,n,n);
% M2 = reshape(m2,n,n);
%res = y - M1*M2*data;
res = res(:);


end