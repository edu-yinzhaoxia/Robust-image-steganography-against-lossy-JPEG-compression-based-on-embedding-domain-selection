
%% 自适应嵌入域

function [cover_round,change_p,change_m,rho_p,rho_m] = I_gmas(cover_Path,rho1_P,rho1_M,tab_m,usable_DCT_num)
    
    cover_spa = imread(cover_Path);
    bits = 8;            %空域像素值转换为非量化DCT系数要将像素值向下平移，使像素值的分布与0对称
    cover_spa = double(cover_spa) - 2^(round(bits) - 1);
    [xm,xn] = size(cover_spa);
    
    % t = dctmtx(8);  % DCT变换矩阵
    % fun = @(xl) (t*xl*(t'));
    % cover_DCT = blkproc(cover_spa,[8 8],fun);  %分块DCT变换
    
    fun = @(x)dct2(x.data);
    cover_DCT = blockproc(cover_spa,[8 8],fun); % dct变化后的系数
    
    m_block = floor(xm/8);
    n_block = floor(xn/8);
    
    %% 广义抖动调制机制计算过程
    G = 1;
    %% usable_DCT_num 初始化嵌入域大小的参数
    
    n_lsb = 0;
    cover_round = zeros(1,m_block*n_block* usable_DCT_num); %
    cover_lsb = zeros(1,m_block*n_block* usable_DCT_num);
    change_p = zeros(1,m_block*n_block* usable_DCT_num); %GMAS 嵌入时非零DCT系数 +- 改变量
    change_m = zeros(1,m_block*n_block* usable_DCT_num);
    rho_p = zeros(1,m_block*n_block* usable_DCT_num);   %GMAS 嵌入时DCT系数 +-1 的失真
    rho_m = zeros(1,m_block*n_block* usable_DCT_num);
    rho1_P_deq  = zeros(1,m_block*n_block* usable_DCT_num); %非量化DCT系数 +-1 的失真
    rho1_M_deq = zeros(1,m_block*n_block* usable_DCT_num);
    
    for bm = 1:m_block
        for bn = 1:n_block
            for i = 1:8
                for j = 1:8
                    %% 这里的循环条件
                    % 循环条件设置
                    if adaptive_condition(i,j,usable_DCT_num)==1
                        n_lsb = n_lsb + 1;
                        rho1_P_deq(n_lsb) = rho1_P( (bm-1)*8+i,(bn-1)*8+j ) / double( tab_m(i,j) ); % 失真/量化表对应的值
                        rho1_M_deq(n_lsb) = rho1_M((bm-1)*8+i,(bn-1)*8+j)/double(tab_m(i,j));
                        yd = cover_DCT((bm-1)*8+i,(bn-1)*8+j);  %非量化的DCT系数
                        tab_q = double(tab_m(i,j))/G;  %量化步长
                        cover_round(n_lsb)=round(yd/tab_q); % 量化和舍入操作后的系数
                        dnum1 = round(yd/tab_q); % 暂存量化舍入系数的奇偶性
                        if mod(dnum1,2)==0
                            cover_lsb(n_lsb)=0; % 偶数，代表0消息
                            dnum2 = floor(yd/tab_q); % 下界
                            if mod(dnum2,2)==1 % 下界是奇数 (说明量化前，系数值，大于0.5)
                                % 若嵌入1，就是左右修改，向右(dnum2+2)，距离就是相减
                                change_p(n_lsb) = (dnum2+2)*tab_q-yd;
                                % 失真修改
                                rho_p(n_lsb) = change_p(n_lsb)*rho1_P_deq(n_lsb);
                                % 向左要到dnum2*tab_q
                                change_m(n_lsb) = -1*(yd-dnum2*tab_q);
                                rho_m(n_lsb) = -1*change_m(n_lsb)*rho1_M_deq(n_lsb);
                            else
                                change_p(n_lsb) = (dnum2+1)*tab_q-yd;
                                rho_p(n_lsb) = change_p(n_lsb)*rho1_P_deq(n_lsb);
                                change_m(n_lsb) = -1*(yd-(dnum2-1)*tab_q);
                                rho_m(n_lsb) = -1*change_m(n_lsb)*rho1_M_deq(n_lsb);
                            end
                        else
                            cover_lsb(n_lsb)=1;
                            dnum2 = floor(yd/tab_q);
                            if mod(dnum2,2)==1
                                change_p(n_lsb) = (dnum2+1)*tab_q-yd;  %加上这么多（加正数）
                                rho_p(n_lsb) = change_p(n_lsb)*rho1_P_deq(n_lsb);
                                change_m(n_lsb) = -1*(yd-(dnum2-1)*tab_q); %加上这么多（加负数）
                                rho_m(n_lsb) = -1*change_m(n_lsb)*rho1_M_deq(n_lsb);
                            else
                                change_p(n_lsb) = (dnum2+2)*tab_q-yd;  %加上这么多（加正数）
                                rho_p(n_lsb) = change_p(n_lsb)*rho1_P_deq(n_lsb);
                                change_m(n_lsb) = -1*(yd-dnum2*tab_q); %加上这么多（加负数）
                                rho_m(n_lsb) = -1*change_m(n_lsb)*rho1_M_deq(n_lsb);
                            end
                        end
                        
                    end
                end % 最内层循环结束 for j
            end
        end
    end
    
end