%% readme
% ��Ҫ�޸�������·����names�������Լ�������ļ���


addpath(fullfile('../JPEG_Toolbox'));
%��ȡn��ͼ�������
%-------------------------------------------------------------------------
tic;
clear;clc;
%% ����
n=2000;%%%
QF=75;
F=zeros(n,8000);
% file_path='F:\codes\data\Imageset\Imageset\BossBase-1.01-cover.tar\jpg75\';% ���壬�����ظ�ʹ��
file_path='../stego_dir/';% ��дͼ��

img_path_list = dir(strcat(file_path,'*.jpg'));%��ȡ���ļ���������jpg��ʽ��ͼ��
% parpool local 12

%% ���߳���
CoreNum=4; %�趨����CPU��������
if isempty(gcp('nocreate')) %�������δ����
    parpool(CoreNum);
end

parfor i=1:n
    image_name = img_path_list(i).name;%ͼ����
    I_STRUCT = jpeg_read(strcat(file_path,image_name));
    F(i,:)=DCTR(I_STRUCT,QF);
    fprintf(['�� ',num2str(i),' ��ͼ��-------- ok','\n']);
end
% parpool close;

load('./names2000.mat');%% ***************
% save(['./dxl_cover_',num2str(n),'_',num2str(QF),'_',datestr(now,29),'.mat'],'F','names');
save(['./dxl_stego_',num2str(n),'_',num2str(QF),'_',datestr(now,29),'.mat'],'F','names');

toc;
%--------------------------------------------------------------------------