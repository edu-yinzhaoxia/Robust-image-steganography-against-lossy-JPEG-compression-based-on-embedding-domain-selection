function result=adaptive_condition(i,j,usable_DCT_num)
    result = 0;
    if usable_DCT_num==64
        %% 全频嵌入，默认直接返回1
        
         result = 1;
        
    end
    
    if usable_DCT_num==35
        if (i+j==3)||(i+j==4)||(i+j==5) ||(i+j==6)||(i+j==7)||(i+j==8)||(i+j==9) %35 个DCT系数
            result = 1;
        end
        
    end
    
    if usable_DCT_num==33
        if (i+j==4)||(i+j==5) ||(i+j==6)||(i+j==7)||(i+j==8)||(i+j==9) 
            result = 1;
        end
        
    end
     
    if usable_DCT_num==30
        if (i+j==5) ||(i+j==6)||(i+j==7)||(i+j==8)||(i+j==9)
            result = 1;
        end
        
    end
        
    if usable_DCT_num==26
        if (i+j==6)||(i+j==7)||(i+j==8)||(i+j==9) 
            result = 1;
        end
        
    end
        
        
    if usable_DCT_num==21
        if (i+j==7)||(i+j==8)||(i+j==9)
            result = 1;
        end
        
    end
    
    
    
    
end