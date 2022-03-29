

load('mass_balance_2D_inlet.mat')
out = output{1};
load('mass_balance_2D')
out1 = output{1};

t = b_struct.dt*0.2*(1:500);
%plot1
d = double(diff(out.x_s_save,1,2));
plot(t(1:499),[d(500,:); max(d) ;mean(d); min(d)]./10'), hold on
plot(t(1:499),diff(mean(out1.x_s_save(:,1:500)))./10), hold on

figure
b = out.x_s_save-1800;
plot(t,[b(500,:); max(b); mean(b) ;min(b)]'), hold on
plot(t,mean(out1.x_s_save(:,1:500))-1800), hold on


figure,
plot(linspace(0,5,100000),out.Qinlet.*b_struct.dt./b_struct.dy,'Color',[0.8 0.8 0.8]), hold on
plot(t,out.Qoverwash(1:200:end)./b_struct.dy), hold on
%plot(t,out.Qinlet(1:200:end))

s = smooth(out.Qinlet,50/b_struct.dt);
plot(t,s(1:200:end)./b_struct.dy)

%plot(t,sum(reshape(out.Qinlet,200,500))./b_struct.dy)
