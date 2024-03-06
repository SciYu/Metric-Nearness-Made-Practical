function [X] = hlwb_projection(X0, D0, niters)
% function [X] = hlwb_projection(X0, D0, niters)
%
% The alternating projection stage iteratively refnes the 
% approximate solution X0 to the optimal X. (see reference)
%
% @param X0      approximate solution
% @param D0      initial non-metric distance matrix
% @param niters  default 100
%
% @return X      optimal solution to metric nearness problem
%
% <Reference>
% [1] Wenye Li, Fangchen Yu, and Zichen Ma. "Metric nearness made practical." AAAI, 2023.

if (nargin < 3)
    niters = 100;
end

n = size(D0, 1);
t = 2;
lambda = 1/t;
X = X0;
for iter = 1 : niters
    X = lambda * D0 + (1-lambda) * X;
    lambda = 0.382/iter;
    for i = 1 : n
        for j = i+1 : n
            for k = 1 : n
                delta = (X(i,j)-X(i,k)-X(j,k)) / 3;
                if (k~=i) && (k~=j) && (delta>0)
                    X(i,j) = X(i,j) - delta; X(j,i) = X(i,j);
                    X(i,k) = X(i,k) + delta; X(k,i) = X(i,k);
                    X(j,k) = X(j,k) + delta; X(k,j) = X(j,k);
                end
            end
        end
    end
    if mod(iter, 10) == 0
        [~, CSR] = ismetric(X);
        fprintf('\niter=%d',iter);
        fprintf(': CSR=%0.2f',CSR);
    end
end

X = (X+X') / 2;
X(1:n+1:n^2) = 0;

end
