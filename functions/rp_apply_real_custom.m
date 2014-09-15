function G = rp_apply_real_custom(X,proj)

D = size(X,2);

V = proj.W*X;

% Add offset b
% Suboptimal!!! Do preprocessing instead...

for i = 1:D
    V(:,i) = V(:,i) + proj.b;
end

G = sqrt(2/D)*cos(V);
end
