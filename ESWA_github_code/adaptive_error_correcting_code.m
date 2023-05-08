%% 自适应纠错码，和生成纠错编码后的信息

function [nn,kk,mm,rs_encoded_msg] = adaptive_error_correcting_code(nn,kk,mm,raw_msg_len,raw_msg)

    zeros_padding_num = ceil(raw_msg_len/kk/mm)*kk*mm - raw_msg_len; %需要补零的个数
    zeros_padding_msg = zeros(1, raw_msg_len + zeros_padding_num);
    zeros_padding_msg(1:raw_msg_len) = raw_msg;
    zeros_padding_msg(raw_msg_len+1 : raw_msg_len + zeros_padding_num) = 0;
    [rs_encoded_msg] = rs_encode_yxz(zeros_padding_msg,nn,kk);
    
end