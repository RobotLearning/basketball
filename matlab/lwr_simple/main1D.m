clear
close all
dbstop if error

h1 = figure; grid on; hold on;
xlabel 'Input X';
%set_fig_position([0.185 0.317 0.596 0.475]);
set_fig_position([0.259 0.00278 0.692 0.997]);


% create training data 1D
[train.input, train.output] = make_training_data();
[train.input, train.output] = randomize_train(train.input, train.output);
axis equal;
testInVec  = min(train.input)-5:1:max(train.input)+15;

height = max(train.output)*0.5;


for testIn = testInVec
    
    % compute weights for the given query point
    %testIn = 2.5;
    diff = bsxfun(@minus, train.input, testIn);
    kwidth = 0.1;
    w = exp( -kwidth*sum(diff.^2, 2) );
    
    norm_w = height*(w./max(w));

    try
       delete(hc.fig2);
    end    

    hc.fig2(1) = plot(train.input, norm_w, 'r.');
    hc.fig2(2) = plot([testIn testIn]', [0 height]');


    figure(h1)
    xlim([testInVec(1) testInVec(end)]); 
    try
        delete(hc.fig1);
    end
    % reweight input and output training data
    [testOutMean, testOutCov] = w_normal(w, train.input, train.output, testIn);
    
    hc.fig1(1) = plot(testIn, testOutMean, 'mo', 'MarkerSize', 15, 'MarkerFaceColor', 'none');
    plot(testIn, testOutMean, 'Color', [.5 .5 .5], 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'none');
    hc.fig1(2) = plot([testIn testIn],...
        [testOutMean-2*sqrt(testOutCov) testOutMean+2*sqrt(testOutCov)], 'm', 'LineWidth', 5, 'MarkerFaceColor', 'none');
    plot([testIn testIn],...
        [testOutMean-2*sqrt(testOutCov) testOutMean+2*sqrt(testOutCov)], 'Color', [.5 .5 .5], 'LineWidth', 1, 'MarkerFaceColor', 'none');
    
    drawnow;
    pause(.15);
    
end
try
    delete(hc.fig1);
end























            