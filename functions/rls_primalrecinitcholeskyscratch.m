function [cfr] = rls_primalrecinitcholeskyscratch(X, y, opt)
% rls_primalrecinitcholeskyscratch(X,y,opt)
% computes a classifier for the primal formulation of RLS using the Cholesky decomposition of the covariance matrix. 
% The variables necessary for further recursive update are stored in the
% output structure.
% The regularization parameter is set to the one found in opt.paramsel.
% In case of multiclass problems, the regularizers need to be combined with the opt.singlelambda function.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- paramsel.lambdas (set by the paramsel_* routines)
%   fields with default values set through the defopt function:
%		- singlelambda
%       - opt.randfeats.D
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -b: parameter vector
% -C: empty matrix
% -X: empty matrix
% -R: Cholesky factor of the covariance matrix C
% -XtX: empty matrix

if isfield(opt.paramsel,'lambdas')
    lambda = opt.singlelambda(opt.paramsel.lambdas);
else
    error('ERROR: No regularization parameter was found in opt.paramsel.lambdas.')
end

% Compute R from starting from XtX
d = size(X,2);
n = size(y,1);


% if isfield(opt,'kernel');
%     if isfield(opt.kernel,'XtX');
%         XtX = opt.kernel.XtX;
%     else
%         XtX = X'*X;
%     end
%     if isfield(opt.kernel,'Xty');
%         Xty = opt.kernel.Xty;
%     else
%         Xty = X'*y;
%     end
%     if isfield(opt.kernel,'n');
%         n = opt.kernel.n;
%     else
%         n = size(y,1);
%     end
% else
%     XtX = X'*X;
%     Xty = X'*y;
% end

% % This instruction computes R from starting from XtX
% R = chol(XtX + (n*lambda)*eye(d));
% 
% % Initialize parameter vector b and weights
% b = Xty;    % From previous training
% %W = rls_primal_driver( XtX, Xty, n, lambda );   % From previous matrix
% W = R\(R'\Xty);   % From previous matrix

% Compute R, W, b from scratch, disregarding X
    b = zeros(d,size(y,2));    % From scratch
    W = zeros(d,size(y,2));    % From scratch
    R = sqrt(n*lambda)*eye(d);

% Save in OPT
cfr.W = W;
cfr.b = b;
cfr.C = [];
cfr.X = [];
cfr.R = R;  % Store R instead of Cinv
cfr.XtX = [];  % DEBUG: Return XtX    

