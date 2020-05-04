function [Pass_Fail_Set] = PIIQ_v3(m,y_int,x_random,y_random)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    Pass_Fail_Set_3=zeros(1,4);
    x_rand=x_random;
    y_rand=y_random;
    for i=1:4
        test_pts=[m(i)*x_rand+y_int(i),(y_rand-y_int(i))/m(i)]; %y_test point, x_test point
        if i==1
            if x_rand<=test_pts(2) && test_pts(1)<=y_rand
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==2
            if x_rand<=test_pts(2) && test_pts(1)>=y_rand
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==3
            if x_rand>=test_pts(2) && test_pts(1)>=y_rand
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
       if i==4
            if x_rand>=test_pts(2) && test_pts(1)<=y_rand
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
       end
       Pass_Fail_Set_3(i)=Pass_Fail;
    end
    Pass_Fail_Set=Pass_Fail_Set_3;
end

