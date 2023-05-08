%% readme
% 需要修改数量，路径，names数量，以及保存的文件名


addpath(fullfile('../JPEG_Toolbox'));
%提取n幅图像的特征
%-------------------------------------------------------------------------
tic;
clear;clc;
%% 参数
n=2000;%%%
QF=75;
F=zeros(n,8000);
% file_path='F:\codes\data\Imageset\Imageset\BossBase-1.01-cover.tar\jpg75\';% 载体，可以重复使用
file_path='../stego_dir/';% 隐写图像

img_path_list = dir(strcat(file_path,'*.jpg'));%获取该文件夹中所有jpg格式的图像
% parpool local 12

%% 多线程跑
CoreNum=4; %设定机器CPU核心数量
if isempty(gcp('nocreate')) %如果并行未开启
    parpool(CoreNum);
end

parfor i=1:n
    image_name = img_path_list(i).name;%图像名
    I_STRUCT = jpeg_read(strcat(file_path,image_name));
    F(i,:)=DCTR(I_STRUCT,QF);
    fprintf(['第 ',num2str(i),' 幅图像-------- ok','\n']);
end
% parpool close;

load('./names2000.mat');%% ***************
% save(['./dxl_cover_',num2str(n),'_',num2str(QF),'_',datestr(now,29),'.mat'],'F','names');
save(['./dxl_stego_',num2str(n),'_',num2str(QF),'_',datestr(now,29),'.mat'],'F','names');

toc;
%--------------------------------------------------------------------------