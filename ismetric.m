function [metric, CSR] = ismetric(D, epsilon)
% function [metric, CSR] = ismetric(D, epsilon)
%
% Check if a given matrix meet the distance metric requirements
%
% @param D         The input distance matrix
% @param epsilon   Preferred precision, default 1e-3 
%
% @return metric   metric = 1 if D is a metric, otherwise metric = 0
% @return CSR      Constraints Satisfaction Ratio (CSR) for triangle inequalities
%
% <Reference>
% [1] Wenye Li, Fangchen Yu, and Zichen Ma. "Metric nearness made practical." AAAI, 2023.

if (nargin < 2)
    epsilon = 1e-4;
end
n = size(D, 1);
metric = true;

%% check diagonal
if max(abs(diag(D))) > epsilon
    metric = false;
    CSR = NaN;
    return
end

%% check symmetry
if ~issymmetric(D)
    metric = false;
    CSR = NaN;
    return
end

%% check triangle inequalities
n_total = nchoosek(n,3) * 3; % total number of triangle inequalities
n_fail = 0;                  % the number of violated triangle inequalities
for i = 1 : n-1
    for j = i+1 : n
        num = nnz(D(i,j)-epsilon > D(:,i)+D(:,j));
        if num > 0
            metric = false;
            n_fail = n_fail + num;
        end
    end
end
CSR = 1 - n_fail / n_total; % CSR = 1 - Ratio of violated triangle inequalities

end