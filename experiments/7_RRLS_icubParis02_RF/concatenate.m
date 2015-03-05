% Separate features from labels

deltaHorizon = 50;

% % Concatenate positions, without f/t features
% XsetConc = zeros(totalSize - deltaHorizon , numJoints * deltaHorizon + 2*numJoints);
% 
% for i = 1:size(XsetConc)
%     for j = 0:(deltaHorizon-1)
%         % Add features to be concatenated
%         for k = 1:numJoints
%             XsetConc(i , j * numJoints + k) = Xset(i + j , k);
%         end
%     end
%     % Add other features
%     XsetConc(i , (end-2*numJoints):end) = Xset(i + deltaHorizon - 1 , (end-2*numJoints):end);
% end

% Concatenation (positions, velocities & accelerations)
XsetConc = zeros(totalSize - deltaHorizon , numJoints * deltaHorizon * 3);
for i = 1:size(XsetConc)
    for j = 0:(deltaHorizon-1)
        % Add features to be concatenated
        for k = 1:3*numJoints
            XsetConc(i , j * 3 * numJoints + k) = Xset(i + j , k);
        end
    end
end

Xset = XsetConc;
XsetConc = [];

yset = yset(deltaHorizon : end , :);

% % Concatenation (just positions)
%XsetConc = zeros(size(Xset,1) - deltaHorizon , numJoints * (deltaHorizon + 2));

% for i = 1:size(XsetConc,1)
%     for j = 1 : deltaHorizon
%         % Fill
%         for k = 0 : (numFeats-1)
%             XsetConc(i , (j-1)*numFeats + 1 + k) = Xset(i + deltaHorizon - j , k+1);
%         end
%     end
% end
% 
% % Concatenation (positions, velocities & accelerations)
% XsetConc = zeros(size(Xset,1) - deltaHorizon , numJoints * deltaHorizon * 3);
% for i = 1:size(XsetConc,1)
%     for j = 1 : deltaHorizon
%         % Fill
%         for k = 0 : (numFeats-1)
%             XsetConc(i , (j-1)*numFeats + 1 + k) = Xset(i + deltaHorizon - j , k+1);
%         end
%     end
% end
% 
% Xset = XsetConc;

% yset = [ zeroes(size(yset),deltaHorizon) , yset(deltaHorizon : end , :) ];
