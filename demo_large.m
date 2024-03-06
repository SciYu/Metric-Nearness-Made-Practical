clear all; clc;
addpath('method');

%% load distance matrices
% Dtrue: ground-truth distance matrix
% Dnoise: noisy non-metric distance matrix
load data/sample_large.mat;
fprintf("Demo on Large Data of AAAI-2023 Paper 'Metric Nearness Made Practical'.")
fprintf("\nData is loaded.");

%% Stage I. Embedding calibration and heuristic improvement
Dcal = embedding_calibration(Dnoise);        % Dnoise: noisy non-metric distance matrix
Dheu = heuristic_improve(Dcal, Dnoise, 1);   % 1 iteration is good enough

%% Stage II. HLWB projection with 100 iterations 
Dhlwb = hlwb_projection(Dheu, Dnoise, 100);  % Dheu: starting matrix, 100 iterations

%% Evaluation Metric
NMSE = norm(Dhlwb-Dnoise, 'fro')^2 / norm(Dnoise, 'fro')^2;     % Normalized Mean Squared Error (NMSE)
RSD = norm(Dhlwb-Dtrue, 'fro')^2 / norm(Dnoise-Dtrue, 'fro')^2; % Relative Squared Deviation (RSD)
[metric_old, CSR_old] = ismetric(Dnoise);                       % Constraints Satisfaction Ratio (CSR)
[metric_new, CSR_new] = ismetric(Dhlwb);

if metric_old == 1
    fprintf('\nDnoise is a distance metric.');
else
    fprintf('\nDnoise is a non-metric with CSR = %0.2f.', CSR_old);
end
if metric_new == 1
    fprintf('\nDhlwb is a distance metric.');
else
    fprintf('\nDhlwb is a non-metric with CSR = %0.2f', CSR_new);
end
fprintf('\nNMSE = %0.2f, RSD = %0.2f', NMSE, RSD);
