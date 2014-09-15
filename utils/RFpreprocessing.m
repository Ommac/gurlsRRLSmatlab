function [ XsetRF proj ] = RFpreprocessing( Xset, numRF )

% [ Xset proj ] = RFpreprocessing( Xset, numRF )
% computes the random features projection of the dataset on a
% numRF-dimensional space.
%
% INPUTS:
% - Xset: n X d pattern matrix to be mapped on the RF space
% - numRF: Number of desired random features over whose space Xset has to
%          be projected
%
% OUTPUT: struct with the following fields:
% - XsetRF: Projected dataset
% - proj:
%       .b: numRF realizations of the b random variable
%       .W: numRF realizations of the W random variable

% Extract realizations of the W and b random vars
[proj.W, proj.b] = rp_projections_custom( size(Xset,2), numRF , 'gaussian');

% Generate dataset projection on the random features space
XsetRF = rp_apply_real_custom( Xset', proj )';    % RF mapping
% size(XsetRF)
% size(Xset)