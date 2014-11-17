
if size(RMSEYesRet,1) > 1
    % Plot RMSEs
    figure
    boxplot(mean(RMSEYesRet,3));
    title('RMSE on test set, RLS with retuning');

    figure
    boxplot(mean(RMSENoRet,3));
    title('RMSE on test set, RLS without retuning');

    % Difference
    figure
    boxplot(mean(RMSENoRet - RMSEYesRet,3));
    title('RMSE on test set, RLS without retuning - RLS with retuning');

    %Plot lambdas
    figure;   
    for i = 1:T
        subplot(2,3,i)
        boxplot(lambdasYesRet(:,:,i));
        set(gca,'YScale','log') 
        title(['Selected lambdas, RLS with retuning, output #' num2str(i)]);
    end
    
    figure;    
    for i = 1:T
        subplot(2,3,i)
        boxplot(lambdasNoRet(:,:,i));
        set(gca,'YScale','log') 
        title(['Selected lambdas, RLS without retuning, output #' num2str(i)]);
    end
    
else
    
    % Plot RMSEs
    figure
    bar(RMSEYesRet);
    title('RMSE on test set, RLS with retuning');

    figure
    bar(RMSENoRet);
    title('RMSE on test set, RLS without retuning');

    % Difference
    figure
    bar(RMSENoRet-RMSEYesRet);
    title('RMSE on test set, RLS without retuning - RLS with retuning');

    %Plot lambdas
    figure
    semilogy(lambdasYesRet);
    title('Selected lambdas, RLS with retuning');

    figure
    semilogy(lambdasNoRet);
    title('Selected lambdas, RLS without retuning');
end