function [Pass_Fail_Set] = PIIQ_0_Function(m,y_int,x_rand,y_rand,x,y)

    Pass_Fail_Set=zeros(4,1);
    Pass_Fail=0;
    for ii=1:4
    test_pts=[(y_rand-y_int(ii))/m(ii), m(ii)*x_rand+y_int(ii)];
        if ii==1 
            if abs(m(ii))<1e-3 && y(1)<=y_rand
                Pass_Fail_Set(ii)=1;
            end
        end
        if ii==2
            if abs(m(ii))>1e3 && x(2)>=x_rand
                Pass_Fail_Set(ii)=1;
            elseif x_rand<=test_pts(1) && y_rand<=test_pts(2)
                Pass_Fail_Set(ii)=1;
            end
        end
        if ii==3
            if abs(m(ii))<1e-3 && y(3)>=y_rand
                Pass_Fail_Set(ii)=1;
            elseif x_rand>=test_pts(1) && y_rand<=test_pts(2)
                Pass_Fail_Set(ii)=1;
            elseif m(ii)<0 && x_rand<=test_pts(1) && y_rand<=test_pts(2)
                Pass_Fail_Set(ii)=1;
            end
        end
        if ii==4
            if abs(m(ii))>1e3 && x(4)<=x_rand 
                Pass_Fail_Set(ii)=1;
            elseif x_rand>=test_pts(1) && y_rand<=test_pts(2)
                Pass_Fail_Set(ii)=1;
            end
        end
    end
end



