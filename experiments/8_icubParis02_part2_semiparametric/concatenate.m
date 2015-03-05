% Separate features from labels

deltaHorizon = 50;

% % Concatenation (positions, velocities & accelerations)
% XsetConc = zeros(totalSize - deltaHorizon , numJoints * deltaHorizon * 3);
% for i = 1:size(XsetConc)
%     for j = 0:(deltaHorizon-1)
%         % Add features to be concatenated
%         for k = 1:3*numJoints
%             XsetConc(i , j * 3 * numJoints + k) = Xset(i + j , k);
%         end
%     end
% end

% % Concatenation (positions, velocities, accelerations, F/T pred)
XsetConc = zeros(totalSize - deltaHorizon , (numJoints * 3 + numLabels ) * deltaHorizon );
for i = 1:size(XsetConc)
    for j = 0:(deltaHorizon-1)
        % Add features to be concatenated
        for k = 1:(numJoints * 3 + numLabels )
            XsetConc(i , j * 3 * numJoints + k) = Xset(i + j , k);
        end
    end
end


Xset = XsetConc;
XsetConc = [];

yset = yset(deltaHorizon : end , :);
f_par = f_par(deltaHorizon : end , :);
t_par = t_par(deltaHorizon : end , :);
