function [C] = embedding_calibration(D, mu)
% function [C] = embedding_calibration(D, mu)
%
% Calibrate a distance matrix by an embedding calibration method. (see references)
%
% @param  D   pairwise distance matrix
% @param  mu  default 0.02
%
% @return C   calibrated distance matrix
%
% <References>
% [1] Wenye Li, Fangchen Yu, and Zichen Ma. "Metric nearness made practical." AAAI, 2023.
% [2] Wenye Li, Fangchen Yu. "Calibrating Distance Metrics Under Uncertainty." ECML, 2022.
% [3] Fangchen Yu, et al. "Highly-Efficient Robinson-Foulds Distance Estimation with Matrix Correction." ECAI, 2023.

if nargin < 2, mu = 0.02; end

gamma = -mu / max(D(:));
low = exp(-mu);
K = nearpsd(exp(gamma.*D), 10, low);
C = log(K) ./ gamma;
m = size(C,1);
C(1:m+1:m*m) = 0;
C = (C+C')/2;

end


function [X, iter] = nearpsd(A, maxits, low, high, d)
% function [X, iter] = nearpsd(A, maxits, low, high, d)
%
% Computes the nearest positive semi-definite matrix 
% for a given square matrix.
%
% @param A        a square matrix to be calibrated
% @param maxits   max num of iters allowed, default 100
% @param low      default 0 
% @param high     default 1
% @param d        values of the diagonal elements, default 1
%
% @return X       nearest psd matrix to A
% @return iter    number of iterations taken

if  ~isequal(A,A'), A = (A + A') / 2; end
if nargin < 5, d = 1; end
if nargin < 4, high = 1; end
if nargin < 3, low = 0; end
if nargin < 2, maxits = 100; end

% threshold for convergence & eigs
tolconv = 1.0e-6;
toleigs = 1.0e-5;

n = size(A,1);

U = zeros(n);
Y = A;

[V, D] = eig(Y);
d = diag(D);

iter = 0;
while 1
    T = Y - U;

    % project onto psd matrices
    [Q, D] = eig(T);
    d = diag(D);
    p = d > toleigs*d(n);
    X = Q(:,p) * D(p,p) * Q(:,p)';

    % update correction
    U = X - T;

    % maximum iteration & convergence test
    iter = iter + 1;
    if iter == maxits
        %fprintf('Max iterations reached. ');
        break; 
    end
    if norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
    	break;
    end
    
    % problem-dependent knowledge added here, e.g.
    Y = X;
    Y(1:n+1:n*n) = d;
    Y(Y<low) = low;
    Y(Y>high) = high;
end

Y(1:n+1:n*n) = d;
Y(Y<low) = low;
Y(Y>high) = high;
%fprintf('Number of iterations taken: %4.0f\n',iter);

end

