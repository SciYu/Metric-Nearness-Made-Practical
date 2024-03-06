function [X] = heuristic_improve(x0, D, niter)
% function [X] = heuristic_improve(x0, D, niter)
%
% Starting from x0,
%       min_X || X - D ||^2 s.t. triangular inequalities
% with an improvement heuristic
%
% <Reference>
% [1] Wenye Li, Fangchen Yu, and Zichen Ma. "Metric nearness made practical." AAAI, 2023.

assert(size(D,1)==size(D,2), 'Input D is not square!');
if (nargin < 3); niter = 1; end

X = x0;
n = size(D, 1);

for iter = 1 : niter
    F = X > D;
    for i = 1 : n
        for j = i+1 : n
            d = D(i,j); X(i,j) = d; X(j,i) = d;
            if F(i,j)
                X(i,j) = max(abs(X(:,i)-X(:,j))); 
            else
                X(i,j) = min(X(:,i)+X(:,j));
            end
            X(j,i) = X(i,j);
        end
    end
end

end