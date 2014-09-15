%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           iCubDyn Dataset loading and preprocessing       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = load('Recursive_matlab_functions/data/icubdyn.dat');

% Separate features from labels
numFeats = 12;
numLabels = 6;
[totalSize, totalDimensions] = size(dataset);
Xset = dataset( : , 1:numFeats );
yset = dataset( : , totalDimensions - numLabels + 1 : totalDimensions);

% Normalize dataset
normX = 1/normest( Xset );
Xset = Xset.*normX;

% Set training set
p_tr = 10000;
Xtr_tot = Xset( 1:p_tr , : );
ytr_tot = yset( 1:p_tr , : );
[ntr_tot,d] = size(Xtr_tot);

% Set arrays containing test/recursive-update set
p_te = 10000;
Xte = Xset( (p_tr + 1):(p_te + p_tr) , : );
yte = yset( (p_tr + 1):(p_te + p_tr) , : );

n0 = 3000; %size of first batch to be used for initialization

Xtr = Xtr_tot(1:n0,:);
ytr = ytr_tot(1:n0,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%