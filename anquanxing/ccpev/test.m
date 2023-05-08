% Perform feature extraction and save extracted features
clc;clear;
dbstop if error;
addpath('./jpeg_toolbox');
addpath('./dct_toolbox');

BatchMergedFeatures('test');
fprintf('-----------\nVerification:\n')
TYPE={'cc-merged','merged'};
for k=1:numel(TYPE)
    % Load features that should have been extracted (reference)
    F1=LoadFeatures(['test/data/expected-' TYPE{k} '.fea']);
    % Load extracted features
    F2=LoadFeatures(['test/data/' TYPE{k} '.fea']);
    % Compare both feature sets
    F1=sortrows(F1,size(F1,2));
    F2=sortrows(F2,size(F2,2));
    fprintf('- number of errors in %s: %i\n',TYPE{k},nnz(abs(F1(:)-F2(:))));
end