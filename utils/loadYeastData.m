%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Yeast dataset loading and preprocessing         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(gurls_root,'demo/data/yeast_data.mat'));
Xtr_tot = Xtr;
ytr_tot = ytr;
[ntr_tot,d] = size(Xtr_tot);



n0 = 100; %size of first batch to be used for initialization

Xtr = Xtr_tot(1:n0,:);
ytr = ytr_tot(1:n0,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%