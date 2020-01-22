function [Pass_Fail_Set_Pent] = PI_Pent_Floor_Function(m,y_int,x_random,y_random, x, y)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

    Pass_Fail_Set_Pent=zeros(1,5);
    x_rand=x_random;
    y_rand=y_random;
        for i=1:5
        test_pts=[(y_rand-y_int(i))/m(i), m(i)*x_rand+y_int(i)]; %x_test point, y_test point
        if i==1
            if x_rand>=x(1)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==2 | i==3
            if x_rand<=test_pts(1) && y_rand>=test_pts(2)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==4 | i==5
            if x_rand<=test_pts(1) && y_rand<=test_pts(2)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        Pass_Fail_Set_Pent(i)=Pass_Fail;
        end
end