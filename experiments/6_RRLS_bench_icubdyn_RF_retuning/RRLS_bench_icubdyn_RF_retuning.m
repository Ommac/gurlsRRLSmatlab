% This demo uses the icubdyn data, with random features projections. 
% Parameter selection and initial RLS estimation is carried out on a
% first subset of the training set.

clearAllButBP;
close all;

numExp = 20; % Number of experiment executions

% Subsets size
q = 5000;   % Training subsets size
p = 5000;   % Test subsets size

nLambda = 1000;  % Number of lambda guesses

% Load dataset
%load icubdyn_ABC;
datasetFileName = 'icubdyn.dat';
dataset = load(datasetFileName);
n = size(dataset,1);

% Separate features from labels
numTotSets = fix(n/(p+q));   % Total number of subsets
idxMax = (p+q)*numTotSets;   % Maximum index of the datasset to be considered
numFeats = 12;
numLabels = 6;
[totalSize, totalDimensions] = size(dataset);
X = dataset( 1:idxMax , 1:numFeats );
y = dataset( 1:idxMax , totalDimensions - numLabels + 1 : totalDimensions);
T = size(y,2);

% Normalize dataset
mx = mean(X);
X = X - mx(ones(size(X,1),1),:); 

% Unit norm + scaling
scaling  = 10;
normX = scaling/normest( X );
X = X.*normX;

% Unit variance
% stdx = std(X);
% X = X ./ stdx(ones(size(X,1),1),:); 

%----------------------------------------
%   Random Features Preprocessing
%----------------------------------------

numRF = 200;
[ X, proj ] = RFpreprocessing( X, numRF ) ;

%----------------------------------------

% Prepare training/testing indexes

% Train
idxTr = zeros(q,numTotSets);
for i = 1:q
    idxTr(i,:) = i:(p+q):(idxMax-p-q+i);
end

% Test
idxTe = zeros(p,numTotSets);
for i = 1:p
    idxTe(i,:) = (i+q):(p+q):(idxMax-p+i+q);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% arrays:
                                                 %
yNoRet = zeros(numExp, numTotSets/2 , p , T);    % NoRet output labels
RMSENoRet = zeros(numExp, numTotSets/2,T);       % NoRet RMSEs on test subsets
lambdasNoRet = zeros(numExp, numTotSets/2,T);    % NoRet chosen lambdas for each subset
yYesRet = zeros(numExp, numTotSets/2 , p , T);   % YesRet test output labels
RMSEYesRet = zeros(numExp, numTotSets/2,T);      % YesRet RMSEs on test subsets
lambdasYesRet = zeros(numExp, numTotSets/2,T);   % YesRet chosen lambdas for each subset

for k = 1:numExp
    
    % RLS with no retuning
    RLSNoRetuning;

    % RLS with retuning
    RLSYesRetuning;

end

% Plot

plots_6

% Save profiler results
profsave(profile('info'),'Recursive_matlab_functions/experiments/6_RRLS_bench_icubdyn_RF_retuning/myprofile_results')

