clear all
clc

x=[3.789, 9.146, 3.789, 0, 0];
y=[-3.788 9.146 9.146 5.358 0]+4;

figure(1)
hold on 
grid on
scatter(x,y,'filled','k')
ylim([0 14])

plot([x(1) x(2)],[y(1) y(2)],'k','Linewidth',2)
plot([x(2) x(3)],[y(2) y(3)],'k','Linewidth',2)
plot([x(3) x(4)],[y(3) y(4)],'k','Linewidth',2)
plot([x(4) x(5)],[y(4) y(5)],'k','Linewidth',2)
plot([x(5) x(1)],[y(5) y(1)],'k','Linewidth',2)


m=[(y(2)-y(1))/(x(2)-x(1)), (y(3)-y(2))/(x(3)-x(2)), (y(4)-y(3))/(x(4)-x(3)), (y(5)-y(4))/(x(5)-x(4)), (y(1)-y(5))/(x(1)-x(5))]; %Bottom, Right Top Left
y_int=[y(1)-m(1)*x(1), y(2)-m(2)*x(2), y(3)-m(3)*x(3), y(4)-m(4)*x(4) , y(5)-m(5)*x(5)];
for ii=1:500
    x_rand=rand*14;
    y_rand=rand*14;
    for i=1:5
        test_pts=[(y_rand-y_int(i))/m(i), m(i)*x_rand+y_int(i)]; %x_test point, y_test point
        if i==1
            if y_rand>=test_pts(2) && x_rand<=test_pts(1)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==2
            if y_rand<=y(2)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==3
            if x_rand>=test_pts(1) && y_rand<=test_pts(2)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==4 
            if x_rand>=x(4)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        if i==5
            if x_rand>=test_pts(1) && y_rand>=test_pts(2)
                Pass_Fail=1;
            else
                Pass_Fail=0;
            end
        end
        Pass_Fail_Set_Pent(i)=Pass_Fail;
    end
    if nnz(Pass_Fail_Set_Pent)==5
        scatter(x_rand,y_rand,'g','filled')
    else
        scatter(x_rand,y_rand,'r','filled')
    end
end