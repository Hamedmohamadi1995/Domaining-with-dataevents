
% clc;
% clear;
% close all
% load DH_5
data = Samples(:,3);
%% First Calculate clusteribility index (CI) and relative clasteribily (RC)
DataVar = var(data);                   % Variance of Data
DataRange = range(data);               % Range of Data
CI = (12.*DataVar)./(DataRange.^2);    % Clusteribility Index
R_CI = CI./ min(CI);                   % Relative Clusteribiity Index
%% Second Calculate z transform of data set

NormStand_data = zscore(data); % convert data into normal standard space
NormStand_DataRange = range(NormStand_data);
Rz_min = min(NormStand_DataRange);
%% Third reweighting the data
test = (Rz_min./NormStand_DataRange).*sqrt(R_CI);
z_3 = NormStand_data.*test;
%%


