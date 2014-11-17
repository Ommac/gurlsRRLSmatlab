if(size(avg_nMSE_rec,4)==1)

    figure
    
    plot(avg_nMSE_rec);
    hold all;
    plot(avg_nMSE_batch);

else
    figure

    m = mean(avg_nMSE_rec,4);
    sd = std(avg_nMSE_rec,0,4);
    f = [ m'+2*sd' , flipdim(m'-2*sd',2)]; 
    fill([1:size(avg_nMSE_rec,1) , size(avg_nMSE_rec,1):-1:1] , f, [7 1 7]/8)
    hold on;
    plot(1:size(avg_nMSE_rec,1) , m , 'b' , 'LineWidth',1);
 
    m2 = mean(avg_nMSE_batch,4);
    sd2 = std(avg_nMSE_batch,0,4);
    f2 = [ m2'+2*sd2' , flipdim(m2'-2*sd2',2)]; 
    fill([1:size(avg_nMSE_rec,1) , size(avg_nMSE_rec,1):-1:1] , f2, [1 7 7]/8)
    hold on;
    plot(1:size(avg_nMSE_rec,1) , m2 , 'r' , 'LineWidth',1);
    
    
    % Log scale
    figure

    m = mean(avg_nMSE_rec,4);
    sd = std(avg_nMSE_rec,0,4);
    f = [ m'+2*sd' , flipdim(m'-2*sd',2)];
    semilogx(1:size(avg_nMSE_rec,1) , m , 'b' , 'LineWidth',1);
    hold on;
    h1 = fill([1:size(avg_nMSE_rec,1) , size(avg_nMSE_rec,1):-1:1] , f, 'yellow', 'FaceAlpha', 0.1);
    semilogx(1:size(avg_nMSE_rec,1) , m , 'b' , 'LineWidth',1);    
    alpha(h1,0.1);

 
    m2 = mean(avg_nMSE_batch,4);
    sd2 = std(avg_nMSE_batch,0,4);
    f2 = [ m2'+2*sd2' , flipdim(m2'-2*sd2',2)]; 
    semilogx(1:size(avg_nMSE_rec,1) , m2 , 'r' , 'LineWidth',1);
    hold on;
    h2 = fill([1:size(avg_nMSE_rec,1) , size(avg_nMSE_rec,1):-1:1] , f2, 'c', 'FaceAlpha', 0.1);
    semilogx(1:size(avg_nMSE_rec,1) , m2 , 'r' , 'LineWidth',1); 
        semilogx(1:size(avg_nMSE_rec,1) , m , 'b' , 'LineWidth',1);    

    alpha(h2,0.1);
  
 
end


    
% if numIter == 1 
% 
%     % Batch errors plot
%     figure('Name', 'Batch, RMSE out 1')
%     bar(RMSE_batch(1,:,1))
%     title('Batch, RMSE out 1')
% 
%     figure('Name', 'Batch, RMSE out 2')
%     bar(RMSE_batch(1,:,2))
%     title('Batch, RMSE out 2')
% 
%     figure('Name', 'Batch, RMSE out 3')
%     bar(RMSE_batch(1,:,3))
%     title('Batch, RMSE out 3')
% 
%     figure('Name', 'Batch, RMSE out 4')
%     bar(RMSE_batch(1,:,4))
%     title('Batch, RMSE out 4')
% 
%     figure('Name', 'Batch, RMSE out 5')
%     bar(RMSE_batch(1,:,5))
%     title('Batch, RMSE out 5')
% 
%     figure('Name', 'Batch, RMSE out 6')
%     bar(RMSE_batch(1,:,6))
%     title('Batch, RMSE out 6')
% 
%     figure('Name', 'Batch, average RMSE')
%     bar(mean(RMSE_batch(1,:,:),3))
%     title('Batch, average RMSE')
%     
%     % Recursive errors plot
%     figure('Name', 'Recursive, RMSE out 1')
%     bar(RMSE_rec(1,:,1))
%     title('Recursive, RMSE out 1')
% 
%     figure('Name', 'Recursive, RMSE out 2')
%     bar(RMSE_rec(1,:,2))
%     title('Recursive, RMSE out 2')
% 
%     figure('Name', 'Recursive, RMSE out 3')
%     bar(RMSE_rec(1,:,3))
%     title('Recursive, RMSE out 3')
% 
%     figure('Name', 'Recursive, RMSE out 4')
%     bar(RMSE_rec(1,:,4))
%     title('Recursive, RMSE out 4')
% 
%     figure('Name', 'Recursive, RMSE out 5')
%     bar(RMSE_rec(1,:,5))
%     title('Recursive, RMSE out 5')
% 
%     figure('Name', 'Recursive, RMSE out 6')
%     bar(RMSE_rec(1,:,6))
%     title('Recursive, RMSE out 6')
% 
%     figure('Name', 'Recursive, average RMSE')
%     bar(mean(RMSE_rec(1,:,:),3))
%     title('Recursive, average RMSE')
% else
%     
%     % Batch errors plot
%     figure('Name', 'Batch, RMSE out 1, 20 runs')
%     boxplot(RMSE_batch(:,:,1))
%     title('Batch, RMSE out 1, 20 runs')
% 
%     figure('Name', 'Batch, RMSE out 2, 20 runs')
%     boxplot(RMSE_batch(:,:,2))
%     title('Batch, RMSE out 2, 20 runs')
% 
%     figure('Name', 'Batch, RMSE out 3, 20 runs')
%     boxplot(RMSE_batch(:,:,3))
%     title('Batch, RMSE out 3, 20 runs')
% 
%     figure('Name', 'Batch, RMSE out 4, 20 runs')
%     boxplot(RMSE_batch(:,:,4))
%     title('Batch, RMSE out 4, 20 runs')
% 
%     figure('Name', 'Batch, RMSE out 5, 20 runs')
%     boxplot(RMSE_batch(:,:,5))
%     title('Batch, RMSE out 5, 20 runs')
% 
%     figure('Name', 'Batch, RMSE out 6, 20 runs')
%     boxplot(RMSE_batch(:,:,6))
%     title('Batch, RMSE out 6, 20 runs')
% 
%     figure('Name', 'Batch, average RMSE, 20 runs')
%     boxplot(avg_RMSE_batch)
%     title('Batch, average RMSE, 20 runs')
% 
%     % Recursive errors plot
%     figure('Name', 'Recursive, RMSE out 1, 20 runs')
%     boxplot(RMSE_rec(:,:,1))
%     title('Recursive, RMSE out 1, 20 runs')
% 
%     figure('Name', 'Recursive, RMSE out 2, 20 runs')
%     boxplot(RMSE_rec(:,:,2))
%     title('Recursive, RMSE out 2, 20 runs')
% 
%     figure('Name', 'Recursive, RMSE out 3, 20 runs')
%     boxplot(RMSE_rec(:,:,3))
%     title('Recursive, RMSE out 3, 20 runs')
% 
%     figure('Name', 'Recursive, RMSE out 4, 20 runs')
%     boxplot(RMSE_rec(:,:,4))
%     title('Recursive, RMSE out 4, 20 runs')
% 
%     figure('Name', 'Recursive, RMSE out 5, 20 runs')
%     boxplot(RMSE_rec(:,:,5))
%     title('Recursive, RMSE out 5, 20 runs')
% 
%     figure('Name', 'Recursive, RMSE out 6, 20 runs')
%     boxplot(RMSE_rec(:,:,6))
%     title('Recursive, RMSE out 6, 20 runs')
% 
%     figure('Name', 'Recursive, average RMSE, 20 runs')
%     boxplot(avg_RMSE_rec)
%     title('Recursive, average RMSE, 20 runs')
% 
% end

% % Errors plotting
% 
% ymax =  max(max([RMSE_batch, RMSE_rec]));
% figure('Name', 'RMSE Batch vs Recursive RLS, 7 steps')
% 
% a1 = subplot(1,2,1);
% bar(RMSE_batch);
% title(a1 , '\bfBatch errors');
% axis(a1, [-Inf Inf 0 ymax*1.1 ]);
% 
% a2 = subplot(1,2,2);
% bar(RMSE_rec);
% title(a2 , '\bfRecursive errors');
% axis(a2, [-Inf Inf 0 ymax*1.1]);

%zmax =  max(max([RMSD_batchRLSRF_C(:,:), RMSD_recursive_chol_C(:,:), RMSD_recursive_sherm_C(:,:)]));

% figure('Name','Errors')
% 
% a = subplot(1,2,1);
% bar3(RMSD_batchRLSRF_C);
% title(a , ['\bfBatch RLS',...
%         char(10), ...
%         'RF, Gaussian kernel'])
% %axis(a, [0 size(RMSD_batchRLSRF_C, 2)+0.5 0 size(RMSD_batchRLSRF_C, 1)+0.5 0 zmax]);
%     
% b = subplot(1,2,2);
% bar3(RMSD_recursive_chol_C);
% title(b,['\bfRecursive RLS',...
%         char(10),...
%         'RF, Gaussian Kernel'...
%         char(10),...
%         'Cholesky update'])
% %axis(b, [0 size(RMSD_batchRLSRF_C, 2)+0.5 0 size(RMSD_batchRLSRF_C, 1)+0.5 0 zmax]);
%   
% % Curve fitting
% figure('Name','Curve Fitting')
% % Output 1
% subplot(3,2,1)
% plot(yte(:,1),'DisplayName','y_C','YDataSource','y_C');
% hold all
% plot(optRec.pred(:,1),'DisplayName','y_C recursive','YDataSource','y_C');
% plot(optBatch.pred(:,1),'DisplayName','y_C Batch RLS RandomFeatures (RBF Ker.)','YDataSource','y_C');
% 
% % Output 2
% subplot(3,2,2)
% plot(yte(:,2),'DisplayName','y_C','YDataSource','y_C');
% hold all
% plot(optRec.pred(:,2),'DisplayName','y_C recursive','YDataSource','y_C');
% plot(optBatch.pred(:,2),'DisplayName','y_C Batch RLS RandomFeatures (RBF Ker.)','YDataSource','y_C');
% 
% % Output 3
% subplot(3,2,3)
% plot(yte(:,3),'DisplayName','y_C','YDataSource','y_C');
% hold all
% plot(optRec.pred(:,3),'DisplayName','y_C recursive','YDataSource','y_C');
% plot(optBatch.pred(:,3),'DisplayName','y_C Batch RLS RandomFeatures (RBF Ker.)','YDataSource','y_C');
% 
% % Output 4
% subplot(3,2,4)
% plot(yte(:,4),'DisplayName','y_C','YDataSource','y_C');
% hold all
% plot(optRec.pred(:,4),'DisplayName','y_C recursive','YDataSource','y_C');
% plot(optBatch.pred(:,4),'DisplayName','y_C Batch RLS RandomFeatures (RBF Ker.)','YDataSource','y_C');
% 
% % Output 5
% subplot(3,2,5)
% plot(yte(:,5),'DisplayName','y_C','YDataSource','y_C');
% hold all
% plot(optRec.pred(:,5),'DisplayName','y_C recursive','YDataSource','y_C');
% plot(optBatch.pred(:,5),'DisplayName','y_C Batch RLS RandomFeatures (RBF Ker.)','YDataSource','y_C');
% 
% % Output 6
% subplot(3,2,6)
% plot(yte(:,6),'DisplayName','y_C','YDataSource','y_C');
% hold all
% plot(optRec.pred(:,6),'DisplayName','y_C recursive','YDataSource','y_C');
% plot(optBatch.pred(:,6),'DisplayName','y_C Batch RLS RandomFeatures (RBF Ker.)','YDataSource','y_C');

%% Miscellanea

% Analysis
% 
% % Autocorrelation
% figure('Name','Autocorrelation')
% for i=1:12
%     subplot(3,4,i)
%     [c_ww,lags] = xcorr(X_A(1:15000,i),15000,'coeff');
%     plot(lags,c_ww)
% end

%   
% % Compute MSE
% 
% %mse = res.perf.rmse.^2;
% mse = sum((res.pred-yte).^2)./size(yte, 1);
% 
% % Compute RMSE
% 
% rmse = sqrt(mse);
% 
% % Compute normalized MSE (nMSE) and average nMSE
% %
% %   nMSE: normalized Mean Squared Error (nMSE), de?ned
% %   as the squared error on the test set divided by the 
% %   variance of the target outputs.
% %
% 
% nmse = mse./var(yte,1);
% 
% % Display errors
% disp(sprintf('\nPrediction accurcay is:\n'))
% 
% disp(sprintf(' RMSE \t'))
% disp(sprintf(' %d\t',rmse))
% 
% disp(sprintf('MSE \t'))
% disp(sprintf(' %d\t',mse))
% 
% disp(sprintf(' nMSE \t'))
% disp(sprintf(' %d\t',nmse))
% 
% disp(sprintf(' Average nMSE \t'))
% disp(sprintf(' %d\t', sum(nmse, 2)/size(nmse, 2)))
% 
% % sumTestPred = sum( res.pred( : , : ) );
% % sumTestLabels = sum( yte( : , : ) );
% % nmse = res.perf.rmse.^2 ./ ( sumTestPred .* sumTestLabels )
% % average_nmse = sum(nmse(:))./size( nmse , 2 )
