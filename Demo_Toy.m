%% =================================================================
% This script runs the TWLRS with proposed generalized inverse operation code
%
% Please make sure your data is in range [0, 1].
%
% Reference: Peilin Yang, Yuning Qiu, Zhenhao Huang, Guoxu Zhou, Qibin Zhao
%            "Efficient and Compact Tensor Wheel Decomposition for Tensor Completion"
%
% Created by Peilin Yang (yangpeilin37@gmail.com)
% Sep. 23, 2023
% Updated by 
% July. 25, 2025

%% =================================================================
clc;
clear;
close all;
addpath(genpath([pwd,'/lib']));
addpath(genpath([pwd,'/data']));


X = imread("House.png");
X = double(X)/255;

%% 
sample_ratio = 0.3;
fprintf('### Performing SR: %4.2f ###\n', sample_ratio);
Y_tensor = X;
clear X
Ndim = ndims(Y_tensor);
Nway = size(Y_tensor);
rand('seed',2);
Omega = find(rand(prod(Nway),1)<sample_ratio);
Y_tensor0 = zeros(Nway);
W =zeros(Nway);
Y_tensor0(Omega) = Y_tensor(Omega);
%%
Algorithms = {'Observed','TW-TC'};
EN_TW_TC = 1;
A_num = length(Algorithms);
Re_tensor = cell(A_num,1);
NumIndexes = 2;  
MatrixTimes = zeros(A_num, 1);
alg = 0;

psnr = [];
ssim = [];
rse = [];
fprintf('###################### Please wait......######################\n')
%% Miss_Data
alg = alg+1;
F = zeros(Nway);
F(Omega) = Y_tensor(Omega);
Re_tensor{alg} = F;
[psnr(alg), ssim(alg), rse(alg)] = quality_ypl2(Y_tensor*255, Re_tensor{alg}*255);
%% Perform Our TWLRS_TC Algorithm
alg = alg+1;
if EN_TW_TC
    fprintf('\n');
    disp(['performing ',Algorithms{alg}, ' ... ']);
    % Initialization of the parameters
    % Please refer to our paper to set the parameters
    opts = [];
    opts.tol   = 1e-8;
    opts.maxit = 500;
    opts.rho   = 3;
     opts.lambda= [2, 0.002];
    opts.Xtrue = Y_tensor;
    opts.R = [2, 7, 3; % R_i
              3, 3, 2]; % L_i
    t_s = tic;
    [Re_tensor{alg}, G, C, Out] = TW_LRG_sparsityGc_TC(Y_tensor0, Omega, opts);
     MatrixTimes(alg,:) = toc(t_s);
    [psnr(alg), ssim(alg), rse(alg)] = quality_ypl2(Y_tensor*255, Re_tensor{alg}*255);
end

