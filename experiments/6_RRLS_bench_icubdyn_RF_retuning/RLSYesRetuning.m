% Preform numTotSets/2 separate batch trainings, recomputing the regularization parameter
% at every iteration. Stores the resulting labels and empirical errors

name = 'RLSYesRetuning';
optYesRet = defopt(name);
optYesRet.nlambda = nLambda;
optYesRet.hoperf = @perf_rmse;
%optYesRet.smallnumber= 10^-16;

for i = 1:(numTotSets/2)
    
    optYesRet.seq = {'split:ho','paramsel:hoprimal','rls:primalrecinit'};   % Replace
    %optYesRet.seq = {'split:ho','paramsel:hoprimal','rls:primal'};

    optYesRet.process{1} = [2,2,2];    % Recursive training on A
    optYesRet = gurls(X( idxTr(:,i) , :), y( idxTr(:,i) , :), optYesRet, 1);
    lambdasYesRet(k,i,:) = optYesRet.paramsel.lambdas;    % Save used lambdas

    yYesRet(k,i,:,:) = pred_primal(X(idxTe(:,i),:), y(idxTe(:,i),:), optYesRet);    % Prediction

    %Compute Root Mean Squared Error (RMSE)
    RMSEYesRet(k,i,:) = sqrt(sum((squeeze(yYesRet(k,i,:,:)) - y(idxTe(:,i),:)).^2))./p;

end