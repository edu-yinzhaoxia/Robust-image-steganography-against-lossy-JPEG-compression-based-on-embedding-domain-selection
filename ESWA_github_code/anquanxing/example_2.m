function example_2()

tic;
clear all;
clc;
%%% <ENSEMBLE SETUP> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

settings.cover = './dxl_cover_2000_75_2023-01-24.mat'; % extracted cover  features
settings.stego = './dxl_stego_2000_75_2023-01-27.mat'; % extracted stego  features

% covariance memory caching turned on => speed-up (may require a lot of
% memory; in case of MEMORY issues, disable this by commenting the
% following line)
settings.keep_cov = 1;

number_of_splits = 10; % number of TRN/TST splits

%%% </ENSEMBLE SETUP> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize final results structure
results = cell(number_of_splits,1);

% loop over all TRN/TST splits
for i=1:number_of_splits
    settings.seed_trntst = i; % PRNG seed for i-th TRN/TST split
    results{i} = ensemble(settings); % launch ENSEMBLE classifier
end

% calculate and report mean error over all TRN/TST splits
results_all = zeros(number_of_splits,1); for i=1:number_of_splits, results_all(i) = results{i}.testing_error; end
fprintf('# -------------------------\nAVERAGE TESTING ERROR OVER %i splits: %.4f\n',number_of_splits,mean(results_all));
toc;