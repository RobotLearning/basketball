function [mean, cov] = w_normal(w, input, output, test)
   
obj.X = input;
obj.Y = output;
obj.weights = w;
obj.x_test = test;


M  = size(obj.X, 1);
Dx = size(obj.X, 2);
Dy = size(obj.Y, 2);

% obj.x_test = x_input;
% 
% if isempty(obj.weights)
%     obj.ComputeWeights(x_input);
% end

n = sum( obj.weights);
mu_x = sum( repmat( obj.weights, 1, Dx) .* obj.X, 1) / n;
mu_y = sum( repmat( obj.weights, 1, Dy) .* obj.Y, 1) / n;

x_diff = obj.X - repmat(mu_x, M, 1);
y_diff = obj.Y - repmat(mu_y, M, 1);

cov_x = x_diff'*( repmat( obj.weights, 1, Dx) .* x_diff) / n;
cov_y = y_diff'*( repmat( obj.weights, 1, Dy) .* y_diff) / n;
cov_xy = x_diff'*( repmat( obj.weights, 1, Dy) .* y_diff) / n;


mean = mu_y' + cov_xy' * pinv(cov_x)*(obj.x_test - mu_x)';
cov  = cov_y - cov_xy' * pinv(cov_x)*cov_xy;

if 0
    N = numel(input(:,1));
    
    mu_x = w.*input/N;
    mu_y = w.*output/N;
    
    x_diff = input - mu_x;
    cov_x  = x_diff'*w*x_diff/N;
    
    y_diff = input - mu_y;
    cov_y  = y_diff'*w.*y_diff/N;
    
    cov_xy = x_diff'*w.*y_diff/N;
    
    tx.mu  = mu_x;
    tx.cov = [];
    
    ty.mu = mu_y;
    
end
    
end

% 
% M = size(obj.X, 1);
% Dx = size(obj.X, 2);
% Dy = size(obj.Y, 2);
% 
% obj.x_test = x_input;
% 
% if isempty(obj.weights)
%     obj.ComputeWeights(x_input);
% end
% 
% n = sum( obj.weights);
% mu_x = sum( repmat( obj.weights, 1, Dx) .* obj.X, 1) / n;
% mu_y = sum( repmat( obj.weights, 1, Dy) .* obj.Y, 1) / n;
% 
% x_diff = obj.X - repmat(mu_x, M, 1);
% y_diff = obj.Y - repmat(mu_y, M, 1);
% 
% cov_x = x_diff'*( repmat( obj.weights, 1, Dx) .* x_diff) / n;
% cov_y = y_diff'*( repmat( obj.weights, 1, Dy) .* y_diff) / n;
% cov_xy = x_diff'*( repmat( obj.weights, 1, Dy) .* y_diff) / n;
% 
% mean = mu_y' + cov_xy' * pinv(cov_x) * (obj.x_test - mu_x)';
% cov = cov_y - cov_xy' * pinv(cov_x)* cov_xy;

            
            
            