"Metric Nearnes Made Practical" -- under review of AAAI'2023

The code was implemented and tested in MATLAB R2021a. 
It can be freely distributed and modified for non-commercial purposes.

Main files:
--demo.m: run this demo to calibrate a 200*200 matrix to meet metric conditions.
--demo_large.m: run this demo to calibrate a 1000*1000 matrix to meet metric conditions.

Method files:
--hlwb.c: the HLWB projection algorithm. also see matlab/hlwb.m
--hlwb.mexa64: MEX compiled hlwb.c on linux x64.
--improve.c: heuristic improvement of an approximate metric. also see matlab/improve.m
--improve.mexa64: MEX compiled improve.c on linux x64.
--matlab/hlwb.m: If hlwb.c can't be compiled, use this one.
--matlab/improve.m: If improve.c can't be compiled, use this one.
--nearlaplacian.m: embedding calibration of a noisy metric.
--nearpsd.m: compute the nearest PSD matrix.

Data files:
--data/sample.mat: a 200*200 noisy distance matrix and the ground-truth, can be used for both HLWB and PAF.
--data/sample_trf.txt: the same 200*200 noisy distance matrix in TRF input format.
--data/sample_large.mat: a 1000*1000 noisy distance matrix and the ground-truth, can be used for both HLWB and PAF.
--data/sample_large_trf.txt: the same 1000*1000 noisy distance matrix in TRF input format.

Link of baselines:
--TRF can be downloaded from: https://optml.mit.edu/work/soft/metricn.html
--PAF can be downloaded from: https://github.com/rsonthal/ProjectAndForget
--DeepNorm can be downloaded from: https://github.com/spitis/deepnorms


Sep. 20th, 2022
