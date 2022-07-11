function [] = displayCamera(R,C,col)
%Display cameras
if(~exist('col','var') || isempty(col))
    % plot the axis with different colors
    % z is red
    % y is blue
    % x is green
    for i = 1:3
        if(i==3)
            c = 'r';
        elseif(i==2)
            c = 'b';
        else
            c = 'g';
        end
        quiver3(C(1),C(2),C(3),R(i,1),R(i,2),R(i,3),'Color',c,'LineWidth',1);
        hold on;
    end
else
    % plot cameras with a unique color
     for i = 1:3
         quiver3(C(1),C(2),C(3),R(i,1),R(i,2),R(i,3),'Color',col,'LineWidth',1);
         hold on;
     end
end

