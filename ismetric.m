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
CSR = 1;

%% check diagonal
if max(abs(diag(D))) > epsilon
    metric = false;
    return
end

%% check symmetry
if ~issymmetric(D)
    metric = false;
    return
end

%% check triangular inequalities
count = 0; % the number of violated triangular inequalities
for i = 1 : n-1
    for j = i+1 : n
        num = nnz(D(i,j)-epsilon > D(:,i)+D(:,j));
        if num > 0
            metric = false;
            count = count + num;
        end
    end
end
CSR = 1 - count / (nchoosek(n,3)*3); % CSR = 1 - Ratio of violated triangle inequalities

end