% Computes the realizations of the W and b random variables
% D: # of RF
% d: # of initial features

function [W,b] = rp_projections_custom(d,D,kernel)
switch kernel
    case 'gaussian'
%         W = sqrt(2/D)*randn(D,d);
%         b = sqrt(2/D)*2*pi*rand(D,1);
        W = randn(D,d);
        b = 2*pi*rand(D,1);
    otherwise
        error('cannot sample from that yet.')
end
end
