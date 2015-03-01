% Separate features from labels

deltaHorizon = 50;


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

% Concatenation (positions, velocities & accelerations)
XsetConc = zeros(size(Xset,1) - deltaHorizon , numJoints * deltaHorizon * 3);
for i = 1:size(XsetConc,1)
    for j = 1 : deltaHorizon
        % Fill
        for k = 0 : (numFeats-1)
            XsetConc(i , (j-1)*numFeats + 1 + k) = Xset(i + deltaHorizon - j , k+1);
        end
    end
end

Xset = XsetConc;

% yset = [ zeroes(size(yset),deltaHorizon) , yset(deltaHorizon : end , :) ];
