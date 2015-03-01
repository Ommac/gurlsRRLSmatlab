%% Experiment description
% 
% In this experiment, the full dataset is used to compare the performance of 
% a batch estimator trained on the first p samples with the ones of a recursive
% estimator initialized on the same first p samples and then updated at
% every step.
% 

% Notes:
% - The generalization performance is computed by means of the nMSE
%   (normalized MSE) on each of the 6 outputs.
% - At each time step t, the recursive estimator computes the predictions
%   for all the outputs from 1 to t. Then, it updates the Cholesky factor
%   using the label of the t-th sample. (see Gijsberts et al., "Real-time
%   model learning using Incremental Sparse Spectrum Gaussian
%   Process Regression", Neural Networks, section 4.2.1, 2013)
% - The recursive estimator uses the Cholesky factorization of the
%   coveriance matrix
% - Both estimators use random features for Gaussian kernel approximation

% Initialization

clearAllButBP;

% Dataset parameters
p = 5000;  % Initialization training set size

% Number of random features
%numRF = [200, 500, 1000];
numRF = 2000;

% Load dataset
datasetFileName = 'processedData.csv';
dataset = load(datasetFileName);

% Find non-zero rows
md = mean(dataset,1);
sd = std(dataset,1);

% Separate features from labels
numJoints = 4;
numFeats = 12;
numLabels = 6;
[totalSize, totalDimensions] = size(dataset);
Xset = dataset( : , [10:9+numFeats/3  32+10:32+9+numFeats/3  32*2+10:32*2+9+numFeats/3] );
yset = dataset( : , totalDimensions - numLabels + 1 : totalDimensions);

% Concatenate latest joint positions
concatenate

% Normalize dataset
mx = mean(Xset);
Xset = Xset - mx(ones(size(Xset,1),1),:); 

% Unit norm + scaling
scaling = 1;
normX = scaling/normest( Xset );
Xset = Xset.*normX;

% Unit variance
% stdx = std(Xset);
% Xset = Xset ./ stdx(ones(size(Xset,1),1),:); 

% Set initial training set
XtrInit = Xset( 1:p , : );
ytrInit = yset( 1:p , : );

% Set arrays containing test/recursive-update set
numUpdates = size(Xset,1) - p;
% numUpdates = totalSize-p;

Xte = Xset( p+1:p+numUpdates , : );
yte = yset( p+1:p+numUpdates, : );

% Compute output variance for each output on the test set
outVar = var(yte);

% Declare error arrays
nMSE_batch = zeros(size(Xte,1), numLabels  );
nMSE_rec = zeros(size(Xte,1), numLabels );
avg_nMSE_batch = zeros(size(Xte,1), 1  );
avg_nMSE_rec = zeros(size(Xte,1), 1 );

    
% Fourier Random Features mapping initialization
XtrInitRF = zeros(size(XtrInit,1));
XteRF = zeros(size(Xte,1));

%% Batch RLS with Random Features (RBF kernel), primal formulation
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
optBatch = gurls(XtrInit, ytrInit, optBatch, 1);

% Batch predictions        
optBatch = gurls(Xte, yte, optBatch, 2);

% Error storage
nSE = zeros( 1 , 6 );        
for i = 1 : size(Xte,1)

    % Update normalized squared error
    nSE = nSE + (yte(i,:) - optBatch.pred(i,:)).^2./outVar;
    % Compute and save nMSE
    nMSE_batch( i, :) = nSE/i;

end


%% Recursive RLS, linear regression, primal formulation

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

optRec = gurls(XtrInitRF, ytrInit, optRec, 1);  % Initialization

% Initial training on the first p samples
%         optRec.rls = rls_primalrecupdatecholesky(XtrInitRF, ytrInit, optRec);   % Recursive update step on B

% Predict (out of sample) & update (rank-1)
nSE = zeros( 1 , 6 );
optRec.newprop('pred', zeros(size(yte)));

for i = 1 : size(Xte,1)

    % Prediction
    optRec.pred(i,:) = pred_primal(XteRF(i,:), yte(i,:), optRec);   % Prediction on C

    % Update normalized squared error
    nSE = nSE + ( yte(i,:) - optRec.pred(i,:) ).^2./outVar;
    % Compute and save nMSE
    nMSE_rec( i, : ) = nSE/i;

    optRec.rls = rls_primalrecupdatecholesky(XteRF(i, :), yte(i, :), optRec);   % Recursive update step on B
end

% Separate  forces from torques
f_rec = optRec.pred(:,1:3);
t_rec = optRec.pred(:,4:6);
f_bat = optBatch.pred(:,1:3);
t_bat = optBatch.pred(:,4:6);
nte = size(Xte,1);
f_gnd = yte(:,1:3);
t_gnd = yte(:,4:6);

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

plots_7
