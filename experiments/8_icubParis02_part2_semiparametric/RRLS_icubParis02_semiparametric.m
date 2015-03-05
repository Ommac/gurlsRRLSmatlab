%% Experiment description
%
% Semiparametric inverse dynamics learning

% Initialization
close all;
clearAllButBP;

%% Load ground truth

datasetGTFileName = 'part2.csv';
datasetGT = importdata(datasetGTFileName, ',', 3);
% datasetFileName = 'icubdyn.dat';
% dataset = dlmread(datasetFileName , ',');

skip = 0;

%% Load parametric predictions (CAD & adaptive)
predFileName = 'results.csv';
pred = load(predFileName);

% Adaptive RBD parameters
% f_par = pred(:,1:3);
% t_par = pred(:,4:6);

% CAD RBD parameters
f_par = pred(:,7:9);
t_par = pred(:,10:12);

pred = [];

%% Load data

% Dataset parameters
p = 5000;  % Initialization training set size

% Number of random features
numRF = 500;

% Load dataset
datasetFileName = 'processedData.csv';
dataset = load(datasetFileName);

% Separate features from labels
numJoints = 4;
numFeats = 12;
numLabels = 6;
[totalSize, totalDimensions] = size(dataset);
Xset = dataset( : , [10:9+numFeats/3  32+10:32+9+numFeats/3  32*2+10:32*2+9+numFeats/3] );

% Add current prediction to features
Xset = [ Xset , f_par , t_par];

% Filtered output
yset = dataset( skip+1:end , totalDimensions - numLabels + 1 : totalDimensions);

% Raw output
% yset = datasetGT.data(skip+1:end,13:18);

% Concatenate latest joint positions
concatenate

yset_delta = yset - [f_par , t_par];

% Normalize dataset
mx = mean(Xset);
Xset = Xset - mx(ones(size(Xset,1),1),:); 

% Unit norm + scaling
scaling = 1;
normX = scaling/normest( Xset );
Xset = Xset.*normX;

% Set initial training set
XtrInit = Xset( 1:p , : );
ytrInit_delta = yset_delta( 1:p , : );

% Set arrays containing test/recursive-update set
numUpdates = size(Xset,1) - p;
% numUpdates = 2000;
% numUpdates = 5000;
% numUpdates = totalSize-p;

Xte = Xset( p+1:p+numUpdates , : );
yte = yset( p+1:p+numUpdates, : );
f_gnd = yte(:,1:3);
t_gnd = yte(:,4:6);
yte_delta = yset_delta( p+1:p+numUpdates, : );

%% Batch RLS with Random Features (RBF kernel), primal formulation
%  CAD prediction refinement

% Batch training on the first p samples
% Testing on the rest

name = 'BatchRLSRandFeats';
optBatch = gurls_defopt(name);
optBatch.seq = { 'split:ho' , 'paramsel:horandfeats' , 'rls:randfeats' , ...
    'pred:randfeats' , 'perf:rmse'};

optBatch.randfeats.D = numRF;    % Set the number of random features
optBatch.nlambda = 200;
optBatch.process{1} = [2,2,2,0,0]; % Batch training on A
optBatch.process{2} = [3,3,3,2,2]; % Batch prediciton on B

optBatch.hoperf = @perf_rmse;   % Set performance measure

% Batch RLS run
optBatch = gurls(XtrInit, ytrInit_delta, optBatch, 1);

% Batch predictions        
optBatch = gurls(Xte, yte_delta, optBatch, 2);


%% Recursive RLS, linear regression, primal formulation
%  CAD prediction refinement

% Random features preprocessing
% Instantiate datasets

% Initial training set
XtrInitRF = rp_apply_real(XtrInit',optBatch.rls.proj)';
XteRF = rp_apply_real(Xte',optBatch.rls.proj)';

% Pipeline definition
name = 'RecursiveRLSCholesky';
optRec = gurls_defopt(name);
optRec.newprop('paramsel.lambdas' , optBatch.paramsel.lambdas);

optRec.seq = { 'rls:primalrecinitcholesky' };
%optRec.seq = { 'split:ho' , 'paramsel:hoprimal' , 'rls:primalrecinitcholesky' };

%optRec.process{1} = [2,2,2]; % Parameter selection & initial training of the recursive estimator
optRec.process{1} = [2]; % Parameter selection & initial training of the recursive estimator

optRec = gurls(XtrInitRF, ytrInit_delta, optRec, 1);  % Initialization

% Initial training on the first p samples
%         optRec.rls = rls_primalrecupdatecholesky(XtrInitRF, ytrInit, optRec);   % Recursive update step on B

% Predict (out of sample) & update (rank-1)
nSE = zeros( 1 , 6 );
optRec.newprop('pred', zeros(size(yte_delta)));

for i = 1 : size(Xte,1)

    % Prediction
    optRec.pred(i,:) = pred_primal(XteRF(i,:), yte_delta(i,:), optRec);   % Prediction on C
    
    % Update
    optRec.rls = rls_primalrecupdatecholesky(XteRF(i, :), yte_delta(i, :), optRec);   % Recursive update step on B
end

%% Separate  forces from torques
if exist('deltaHorizon') == 0
    deltaHorizon = 0;
end
% f_rec = optRec.pred(:,1:3) + f_par(p+deltaHorizon+1:p+deltaHorizon+numUpdates,:);
% t_rec = optRec.pred(:,4:6) + t_par(p+deltaHorizon+1:p+deltaHorizon+numUpdates,:);
% f_bat = optBatch.pred(:,1:3) + f_par(p+deltaHorizon+1:p+deltaHorizon+numUpdates,:);
% t_bat = optBatch.pred(:,4:6) + t_par(p+deltaHorizon+1:p+deltaHorizon+numUpdates,:);
f_rec = optRec.pred(:,1:3) + f_par(p+1:p+numUpdates,:);
t_rec = optRec.pred(:,4:6) + t_par(p+1:p+numUpdates,:);
f_bat = optBatch.pred(:,1:3) + f_par(p+1:p+numUpdates,:);
t_bat = optBatch.pred(:,4:6) + t_par(p+1:p+numUpdates,:);

nte = size(Xte,1);


% Compute average RMSE on all the outputs

avg_mse_f_rec = zeros(nte,1);
avg_mse_t_rec = zeros(nte,1);
avg_mse_f_bat = zeros(nte,1);
avg_mse_t_bat = zeros(nte,1);

% Recursive
for i = 1:nte
    if i > 1
        avg_mse_f_rec(i) = (i-1)/i * avg_mse_f_rec(i-1) + mean((f_rec(i,:) - f_gnd(i,:)).^2 / i , 2);
        avg_mse_t_rec(i) = (i-1)/i * avg_mse_t_rec(i-1) + mean((t_rec(i,:) - t_gnd(i,:)).^2 / i , 2);
    else    
        avg_mse_f_rec(i) = mean((f_rec(i,:) - f_gnd(i,:)).^2 , 2);
        avg_mse_t_rec(i) = mean((t_rec(i,:) - t_gnd(i,:)).^2 , 2);        
    end    
end

% Batch
for i = 1:nte
    if i > 1
        avg_mse_f_bat(i) = (i-1)/i * avg_mse_f_bat(i-1) + mean((f_bat(i,:) - f_gnd(i,:)).^2 / i , 2);
        avg_mse_t_bat(i) = (i-1)/i * avg_mse_t_bat(i-1) + mean((t_bat(i,:) - t_gnd(i,:)).^2 / i , 2);
    else    
        avg_mse_f_bat(i) = mean((f_bat(i,:) - f_gnd(i,:)).^2 , 2);
        avg_mse_t_bat(i) = mean((t_bat(i,:) - t_gnd(i,:)).^2 , 2);        
    end    
end

avg_rmse_f_rec = sqrt(avg_mse_f_rec);
avg_rmse_t_rec = sqrt(avg_mse_t_rec);
avg_rmse_f_bat = sqrt(avg_mse_f_bat);
avg_rmse_t_bat = sqrt(avg_mse_t_bat);

% Save profiler results
% profsave(profile('info'),'Recursive_matlab_functions/experiments/5_RRLS_bench_icubdyn_RF_nMSE/myprofile_results')


%% Plots

plots_8
