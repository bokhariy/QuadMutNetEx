clear;

% load dataset
load test_data/GBM_nets.mat 
% mySparseMatrix is the binary matrix of mutations (patients in rows, genes in columns)
% GenesNames is a cell matrix with gene names
% nets is a matrix of one or multiple binary networks added together.
% selectionCDFn is a non-uniform cumulative discribution function for random selection of new genes to be evaluated (selectionCDFu is uniform CDF)

% set parameters
C=2.5; % select larger value of C for larger (more genes) datasets
k=1; 
np=0.3;

% run QuaDMutNetEx
nets = GBM_Net_hint_hi2012 + GBM_Net_irefindex + GBM_Net_multinet;
[objective_function,selectedGenes, namesOfSelectedGenes, subNetworkSolution] = QuadMutNetEx(mySparseMatrix, GenesNames, nets, 10000, 30, C, k, np, selectionCDFn)

% calculate solution quality metrics
solutionMetrics=QuadMutNetExMetricsStruct(selectedGenes,mySparseMatrix,subNetworkSolution,C,k,np)

% EXPECTED RESULTS: 
% note that QuaDMutNetEx is a randomized algorithm,
% so the results may not be identical
%
% Elapsed time is 38.353965 seconds.
% 
% objective_function =
% 
%    37.0000
% 
% 
% selectedGenes =
% 
%            1           7          50         239         824        1067
% 
% 
% namesOfSelectedGenes =
% 
%   1×6 cell array
% 
%     'CDKN2A'    'TP53'    'MDM2'    'MDM4'    'MAPK9'    'RPL11'
% 
% 
% subNetworkSolution =
% 
%    (2,1)        3
%    (3,1)        3
%    (4,1)        2
%    (6,1)        2
%    (1,2)        3
%    (3,2)        3
%    (4,2)        3
%    (5,2)        3
%    (6,2)        3
%    (1,3)        3
%    (2,3)        3
%    (4,3)        3
%    (5,3)        1
%    (6,3)        3
%    (1,4)        2
%    (2,4)        3
%    (3,4)        3
%    (5,4)        1
%    (2,5)        3
%    (3,5)        1
%    (4,5)        1
%    (1,6)        2
%    (2,6)        3
%    (3,6)        3
% 
% 
% solutionMetrics = 
% 
%   struct with fields:
% 
%        qObj: 37
%        lObj: 79
%      covTot: 97
%      covOvr: 18
%        covR: 0.8151
%     excessR: 0.1856
%        mCnt: 115
%        sCnt: 119
%        gCnt: 6