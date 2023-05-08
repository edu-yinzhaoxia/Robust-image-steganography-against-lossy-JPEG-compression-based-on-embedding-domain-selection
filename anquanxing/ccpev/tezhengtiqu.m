
tic;
clear;clc;
addpath(fullfile('../../JPEG_Toolbox'));

%%
QF=65;
% file_path='F:\codes\data\Imageset\Imageset\BossBase-1.01-cover.tar\jpg65\';% 载体，可以重复使用
file_path='../../stego_dir/';% 隐写图像

BatchMergedFeatures(file_path);
F=LoadFeatures('./data/merged.fea');

load('./names10000.mat');
% save(['./dxl_cover_ccpev_',num2str(QF),'_',datestr(now,29),'.mat'],'F','names');
save(['./dxl_stego_ccpev_',num2str(QF),'_',datestr(now,29),'_0.10.mat'],'F','names');

toc;
%--------------------------------------------------------------------------