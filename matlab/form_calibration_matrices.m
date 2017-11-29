%% Form two matrices for estimation

function [M,D] = form_calibration_matrices(x)

M1 = reshape(x(1:6),2,3);
M2 = reshape(x(7:12),2,3);
o1 = x(13:14);
o2 = x(15:16);
theta = x(17:19); % Euler angles
d = x(20:22); % from camera to origin
M = [M1, o1; M2, o2; zeros(1,3), 1];
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