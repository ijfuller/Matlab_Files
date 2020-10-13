function [Pass_Fail_Set] = PIIQ_67_5_Function(m,y_int,x_random,y_random,x,y,El_Nod_Conn_Cover,Global_Nodal_Coor_Cover,jj)
    x_rand=x_random;
    y_rand=y_random;    
    Pass_Fail_Set=zeros(4,1);
    Pass_Fail=0;
    
    for ii=1:4
        test_pts=[(y_rand-y_int(ii))/m(ii), m(ii)*x_rand+y_int(ii)];
            if ii==1 
                if x_rand<=test_pts(1) && y_rand>=test_pts(2) 
                    Pass_Fail_Set(ii)=1;
                end
            end
            if ii==2
                if abs(m(ii))<=1e-3 && y(2) >= y_rand
                    Pass_Fail_Set(ii)=1;
                elseif  x_rand<=test_pts(1) && y_rand<=test_pts(2)
                    Pass_Fail_Set(ii)=1;
                end
            end
            if ii==3
                if abs(m(ii))>=1e-3 && x(3) <= x_rand
                    Pass_Fail_Set(ii)=1;
                elseif x_rand>=test_pts(1) && y_rand<=test_pts(2)
                    Pass_Fail_Set(ii)=1;
                end
            end
            if ii==4
                if abs(m(4))>1e3 && x(4) <= x_rand 
                    Pass_Fail_Set(ii)=1;
                elseif x_rand>=test_pts(1) && y_rand>=test_pts(2)
                    Pass_Fail_Set(ii)=1;
                end
            end
    end
end

