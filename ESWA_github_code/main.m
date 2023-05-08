
clear;
clc;
dbstop if error;
tic;
addpath(fullfile('./JPEG_Toolbox'));
addpath(fullfile('./STC3'));
%% 参数
cover_num = 1000;
cover_QF = 75;
attack_QF = 75;
% facebook=71, wechat=80, wb = 60;
% realChannel_QF = 60;

payload = 0.1;
nn = 31; kk = 29; mm = 5;
usable_DCT_num = 64; % 默认值
usable_DCT_nums=[64,35,33,30,26,21];
error_rate_threshold = 0.0001;
kk_threshold = 7;
start_pic = 1;

%% 图像
cover_dir = 'F:\codes\data\Imageset\Imageset\BossBase-1.01-cover.tar\jpg75'; %载体图像所在文件夹
stego_dir = '.\stego_dir'; if ~exist(stego_dir,'dir'); mkdir(stego_dir); end  %载密图像所在文件夹
afterchannel_stego_dir = '.\afterchannel_stego_dir'; if ~exist(afterchannel_stego_dir,'dir'); mkdir(afterchannel_stego_dir); end  %信道处理后载密图像所在文件夹

%信道处理后载密图像所在文件夹
realChannel_stego_dir= '.\realChannel_stego_dir'; if ~exist(realChannel_stego_dir,'dir'); mkdir(realChannel_stego_dir); end
real_bit_error_rate = ones(1,cover_num); %记录测试图像的误码率
%% 批量化
kk_vec = zeros(1,cover_num);
usable_DCT_num_vec = zeros(1,cover_num);
bit_error_rate = ones(1,cover_num); %记录测试图像的误码率
for i_img = start_pic:cover_num
    
    %% 载体图像
    cover_Path = fullfile([cover_dir,'\',num2str(i_img),'.jpg']);
    stego_Path = fullfile([stego_dir,'\',num2str(i_img),'.jpg']);
    afterchannel_stego_Path = fullfile([afterchannel_stego_dir,'\',num2str(i_img),'.jpg']);
    
    % 模拟真实信道使用
    realChannel_stego_Path = fullfile([realChannel_stego_dir,'\',num2str(i_img),'.jpg']);
    
    C_STRUCT = jpeg_read(cover_Path);
    C_COEFFS = C_STRUCT.coef_arrays{1};
    C_QUANT = C_STRUCT.quant_tables{1};
    
    %% 非对称失真
    [rho1_P, rho1_M] = J_UNIWARD_Asy_cost(cover_Path);
    %% 生成秘密信息
    nzAC = nnz(C_COEFFS) - nnz(C_COEFFS(1:8:end,1:8:end));
    raw_msg_len = ceil(payload*nzAC);
    raw_msg = round( rand(1,raw_msg_len) );
    
    %% *************设计循环
    %% 溢出处理
    overflowFlag = 0;
    for i=1:6 % 嵌入域循环
        %         有好的结果记得停止循环
        usable_DCT_num=usable_DCT_nums(i);
        % kk 重置,有条件的重置，否则结束运行，结果丢失kk
        %% 自适应嵌入域
        [cover_round, change_p, change_m, rho_p, rho_m] = I_gmas(cover_Path, rho1_P, rho1_M, C_QUANT,usable_DCT_num);
        
        kk = 29;
        while(bit_error_rate(1,i_img) > error_rate_threshold && kk >= kk_threshold)
            %% 自适应纠错码
            [nn,kk,mm,rs_encoded_msg]=adaptive_error_correcting_code(nn,kk,mm,raw_msg_len,raw_msg);
            
            
            %% STC嵌入
            [suc, stc_n_msg_bits] = stc3_embed_all(rs_encoded_msg, cover_Path, cover_round, rho_p, rho_m, change_p, change_m, cover_QF, stego_Path,usable_DCT_num);
            
            %% 模拟压缩
            imwrite(imread(stego_Path),afterchannel_stego_Path,'quality',attack_QF);
            %% 提取
            [stc_decoded_msg] = stc3_extract_all(afterchannel_stego_Path, stc_n_msg_bits, C_QUANT,usable_DCT_num);
            
            [rs_decoded_msg] = rs_decode_yxz(double(stc_decoded_msg), nn, kk);
            extract_raw_msg = rs_decoded_msg(1:raw_msg_len); %去掉补零
            
            
            if numel(rs_encoded_msg)~=numel(stc_decoded_msg) && overflowFlag==0 % 没溢出前
                overflowFlag =1; % 出现溢出
                kk = 15;% 最小纠错能力
                usable_DCT_num = 21; %最大嵌入域,结果显示鲁棒性降低比较明显，还是约束嵌入域吧。
                continue;
            end
            
            
            %% 计算错误率
            bit_error = double(raw_msg) - double(extract_raw_msg);
            bit_error_number = sum(abs(bit_error));
            bit_error_rate(1,i_img) = bit_error_number/raw_msg_len;
            
            
            if overflowFlag==1 % 前提是信息长度一致，只要信息长度不一致，要么缩小嵌入域（第一次），要么报错（第二次）。
                disp("内循环停止")
                break; % 得出结果就停止循环
            end
            
            % 为了让kk记录纠错能力
            if bit_error_rate(1,i_img) > error_rate_threshold
                kk=kk-2;% 循环条件
            end
            
            
        end
        
        if overflowFlag==1 % 出现问题的循环，能结束就结束
            disp("外循环停止")
            break;
        end
        
        
        if bit_error_rate(1,i_img) <= error_rate_threshold
            break;
        end
        
    end
    % 记录嵌入方法
    kk_vec(1,i_img)=kk;
    usable_DCT_num_vec(1,i_img)=usable_DCT_num;
    
    
    % facebook， WeChat
    % 这里打断点暂停，换入压缩后的图像
    %         微信朋友圈灰度图像，三个通道是相同的信息，截取一个就好
    % 微博
    %         imwrite(imread(stego_Path),realChannel_stego_Path,'quality',realChannel_QF);
    
    %
    %     [stc_decoded_msg] = stc3_extract_all(realChannel_stego_Path, stc_n_msg_bits, C_QUANT,usable_DCT_num);
    %     [rs_decoded_msg] = rs_decode_yxz(double(stc_decoded_msg), nn, kk);
    %     extract_raw_msg = rs_decoded_msg(1:raw_msg_len); %去掉补零
    
    %% 计算真实信道错误率
    %     bit_error = double(raw_msg) - double(extract_raw_msg);
    %     bit_error_number = sum(abs(bit_error));
    %     real_bit_error_rate(1,i_img) = bit_error_number/raw_msg_len;
    %     fprintf('%s\n',['payload: ',num2str(payload),'  image_number: ',num2str(i_img),'  real_bit_error_rate: ',num2str(real_bit_error_rate(1,i_img))]);
    
    
    %  输出每张图像的误码率
    fprintf('%s\n',['payload: ',num2str(payload),'  image_number: ',num2str(i_img),'  error_rate: ',num2str(bit_error_rate(1,i_img))]);
    
end
%  输出所有图像的平均误码率
ave_error_rate = mean(bit_error_rate);
% ave_real_error_rate = mean(real_bit_error_rate);
fprintf('%s\n',['payload: ',num2str(payload),'  ave_error_rate: ',num2str(ave_error_rate)]);
toc;
% fprintf('%s\n',['payload: ',num2str(payload),'  ave_real_error_rate: ',num2str(ave_real_error_rate)]);
% save(['./result_mat/',num2str(cover_num),'_',num2str(payload),'_',num2str(cover_QF),'.mat'],'kk_vec','usable_DCT_num_vec','bit_error_rate');

