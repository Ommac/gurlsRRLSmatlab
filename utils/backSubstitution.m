function w = backSubstitution(R,a,d,T)
% Solving an upper triangular system by back-substitution
% Input matrix U is an n by n upper triangular matrix
% Input vector b is n by 1
% Input scalar n specifies the dimensions of the arrays
% Output vector x is the solution to the linear system
% U x = b
% K. Ming Leung, 01/26/03

w = zeros(n,T);
%for k = 1:T
%     for j = n:-1:1
%         if (U(j,j) == 0)
%             error('Matrix is singular!'); end;
% 
%         w(j,:) = a(j,:)/R(j,j);
%         
%         a(1:j-1,:) = a(1:j-1,:) - R(1:j-1,j)*w(j,:);
%     end
%end

for j = n:-1:1
    if (U(j,j) == 0)
        error('Matrix is singular!'); end;

    
    a(1:j-1,:) = a(1:j-1,:) - R(1:j-1,j)*w(j,:);
    
    w(j,:) = a(j,:)/R(j,j);

end