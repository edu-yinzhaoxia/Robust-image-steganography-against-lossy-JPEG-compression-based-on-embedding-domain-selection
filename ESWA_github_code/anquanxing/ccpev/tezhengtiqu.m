
tic;
clear;clc;
addpath(fullfile('../../JPEG_Toolbox'));

%%
QF=65;
% file_path='F:\codes\data\Imageset\Imageset\BossBase-1.01-cover.tar\jpg65\';% ���壬�����ظ�ʹ��
file_path='../../stego_dir/';% ��дͼ��

BatchMergedFeatures(file_path);
F=LoadFeatures('./data/merged.fea');

load('./names10000.mat');
% save(['./dxl_cover_ccpev_',num2str(QF),'_',datestr(now,29),'.mat'],'F','names');
save(['./dxl_stego_ccpev_',num2str(QF),'_',datestr(now,29),'_0.10.mat'],'F','names');

toc;
%--------------------------------------------------------------------------