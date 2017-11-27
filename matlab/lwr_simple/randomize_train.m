function [input, output] = randomize_train(input, output)

    index  = randperm(length(input(:,1)));
    input  = input(index,:);
    output = output(index,:);

end