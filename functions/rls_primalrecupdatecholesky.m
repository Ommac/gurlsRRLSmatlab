function [rls] = rls_primalrecupdatecholesky (X, y, opt)

% rls_primalrecupdatecholesky(X,y,opt)
% computes a classifier for the primal formulation of RLS, using a
% recursive update of the right Cholesky factor R of matrix C,
% starting from an initial R found in OPT.RLS.
%
% INPUTS:
% -X: input data matrix
% -y: labels matrix
% -OPT: struct of options with the following fields:
%   fields that need to be set through previous gurls tasks:
%		- rls.W (set by rls_primalrecinitcholesky)
%       - rls.R (set by rls_primalrecinitcholesky)
%       - rls.b (set by rls_primalrecinitcholesky)
% 
%   For more information on standard OPT fields
%   see also defopt
% 
% OUTPUT: struct with the following fields:
% -W: matrix of coefficient vectors of rls estimator for each class
% -Cinv: inverse of the regularized kernel matrix in the primal space
% -C: empty matrix
% -X: empty matrix

W = opt.rls.W;
R = opt.rls.R;
b = opt.rls.b;

n = size(X,1);
%d = size(X,2);
%T = size(y,2);

% Sequence of rank-1 updates by application of Cholesky rank-1 updates
for i = 1:n;
    
    % Update b
    b = b + X(i,:)'*y(i,:);
    
    % Update Cholesky factor R
    R = cholupdate(R,X(i,:)');
    
    % TEMP: rewrite for efficiency.
    %W = mldivide(R,(mldivide(R',b)));
end

% Using MATLAB's backslash
% tic
% a = R'\b;
% toc
% tic
% W = R\a;
% toc

% Short form
% tic
W = R\(R'\b);
% toc

% TEST: Directly call backsubstitution function (no tests on R for algorithm
% choice)

% tic
% a = backSubstitution(R',b,d,T);
% toc
% tic
% W2 = backSubstitution(R,a,d,T);
% toc
% 
% W-W2

rls.b = b;
rls.W = W;
rls.C = [];
rls.X = [];
rls.XtX = [];
rls.R = R;

