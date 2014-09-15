% Preform numTotSets/2 separate batch trainings, using the regularization
% parameter chosen at the 1st iteration. Stores the resulting labels and
% empirical errors

name = 'RLSNoRetuning';
optNoRet = defopt(name);
optNoRet.nlambda = nLambda;
optNoRet.hoperf = @perf_rmse;
%optNoRet.smallnumber= 10^-16;
firstRun = 1;

for i = 1:(numTotSets/2)
    if (firstRun == 1)
        optNoRet.seq = {'split:ho','paramsel:hoprimal','rls:primalrecinit'};% Replace
        %optNoRet.seq = {'split:ho','paramsel:hoprimal','rls:primal'};
        optNoRet.process{1} = [2,2,2];    % Recursive training on A
        firstRun = 0;
    else
        optNoRet.paramsel.lambdas = optNoRet.paramsel.lambdas;  % Keep the initial lambda!
        optNoRet.seq = {'rls:primalrecinit'};% Replace
        %optNoRet.seq = {'rls:primal'};
        optNoRet.process{1} = 2;    % Recursive training on A
    end

    optNoRet = gurls(X( idxTr(:,i) , :), y( idxTr(:,i) , :), optNoRet, 1);
    lambdasNoRet(k,i,:) = optNoRet.paramsel.lambdas;    % Save used lambdas
    yNoRet(k,i,:,:) = pred_primal(X(idxTe(:,i),:), y(idxTe(:,i),:), optNoRet);    % Prediction

    %Compute Root Mean Squared Error (RMSE)
    RMSENoRet(k,i,:) = sqrt(sum((squeeze(yNoRet(k,i,:,:)) - y(idxTe(:,i),:)).^2))./p;

end