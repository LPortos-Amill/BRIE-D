function [x_all,y_all] = colon_all(x1,x2,y1,step)

x_all = zeros(length(x1),max(x2-x1)+1);
y_all = zeros(length(x1),max(x2-x1)+1);
for ii=1:length(x1),
    if x2(ii)<x1(ii), continue, end
    x_all(ii,1:(x2(ii)-x1(ii)+1)) = x1(ii):x2(ii);
    y_all(ii,1:(x2(ii)-x1(ii)+1)) = ones(1,(x2(ii)-x1(ii)+1)).*y1(ii);

end
y_all = y_all(:); y_all(x_all==0) = [];
x_all = x_all(:); x_all(x_all==0) = [];

end

