%% Function that is used to fit matrices in calibration

% x are the parameters
% y is the camera output
% residual is the vector of differences between y and predicted output
function res = fit_calibration_matrix(y,data,x)

M = reshape(x,3,5);
res = y - M*data;
% m1 = x(1:n*n);
% m2 = x(n*n+1:end);
% M1 = reshape(m1,n,n);
% M2 = reshape(m2,n,n);
%res = y - M1*M2*data;
res = res(:);


end