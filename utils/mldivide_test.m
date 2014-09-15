range = 1:2000;
time = zeros(1,size(range,2));
for  i = 1:size(range,2)
    
    sz = range(i);
    X = rand(sz);
    U = triu(X);
    b = zeros(sz,1);

    tic;
    W = U\b;    
    time(i) = toc;
end

fitpoly2 = fit((1:size(range,2))',time','poly2');
fitpoly3 = fit((1:size(range,2))',time','poly3');

% Plot the fit
figure
plot(fitpoly2,(1:size(range,2))',time')
fitpoly2

figure
plot(fitpoly3,(1:size(range,2))',time')
fitpoly3