%% barrier volume through time
clear all
load mass_balance_inlet


for ii=1:length(param1),
    

    for jj=1:length(param2),
        out = output{ii,jj};
        Vbarrier = b_struct.dy.*(sum(out.d_sf.*out.d_sf./out.s_sf_save./2.*(1-(b_struct.s_background./out.s_sf_save)),1)+sum(out.h_b_save.*double(out.x_b_save-out.x_s_save),1));
        
        blub(ii,jj) = Vbarrier(end);
       
        Qinlet(ii,jj) = mean(out.Qinlet);
        Qoverwash(ii,jj) = mean(out.Qoverwash);
 %      plot(out.Qoverwash)
    end
end

Qoverwash(:,1) = Qoverwash(:,1)*3; %correction for doing part of the run
F = Qinlet./(Qinlet+Qoverwash);

%figure
%histogram(F(:),[0:0.02:1])
%
% F as a function of grid/timestep
c = {Qinlet./Qinlet(6,4),Qoverwash./Qoverwash(6,4),F./F(6,4),bsxfun(@rdivide,blub,mean(blub))};
for ii=1:4,
figure, hold on
imagesc(c{ii},[0.8 1.2])
colormap(jet)
[a,b] = meshgrid(1:length(param1),1:length(param2));
text(a(:),b(:),cellfun(@(x) num2str(x,3),(num2cell(c{ii}(:))),'UniformOutput',false))
set(gca,'XTick',1:length(param1),'YTick',1:length(param2))
set(gca,'XTickLabel',cellstr(num2str(param1')),'YTicklabel',cellstr(num2str(param2')))
end

%% comparison with same run
load('mass_balance_duplicate.mat')
hold on
clear Qinlet Qoverwash
for ii=1:length(param1),
    for jj=1:length(param2),
        out = output{ii,jj};
        Vbarrier = b_struct.dy.*(sum(out.d_sf.*out.d_sf./out.s_sf_save./2.*(1-(b_struct.s_background./out.s_sf_save)),1)+sum(out.h_b_save.*double(out.x_b_save-out.x_s_save),1));
        Qinlet(ii,jj) = mean(out.Qinlet);
        Qoverwash(ii,jj) = mean(out.Qoverwash);
    end
end

F = Qinlet./(Qinlet+Qoverwash);
histogram(F(:),[0:0.02:1])
%% comparison with just barrier model

q = {'mass_balance_2D_inlet','mass_balance_2D','mass_balance_1D','mass_balance_2D_inlet_no_ff','mass_balance_2D_inlet_with_ff'};
for kk=1:5
load(q{kk})
if kk==3, jj=1; ii=4; else, jj=1; ii=1; end
hold on

        out = output{ii,jj};
        Vbarrier = b_struct.dy.*(sum(out.d_sf.*out.d_sf./out.s_sf_save./2.*(1-(b_struct.s_background./out.s_sf_save)),1)+sum(out.h_b_save.*double(out.x_b_save-out.x_s_save),1));
        
       
        
        plot(Vbarrier./mean(Vbarrier))

        if kk==1, plot(smooth(out.Qinlet(1:end-1),1000)), end
        %Qoverwash(ii,jj) = mean(out.Qoverwash)./param1(ii);
        %plot(out.Qoverwash)

%{
figure
pcolor(Qoverwash)
colormap(jet)
[a,b] = meshgrid(1:5,1:6);
text(a(:),b(:),cellfun(@(x) num2str(x,2),(num2cell(c{ii}(:)./c{ii}(2,2))),'UniformOutput',false))
set(gca,'XTick',1:length(param1),'YTick',1:length(param2))
set(gca,'XTickLabel',cellstr(num2str(param1')),'YTicklabel',cellstr(num2str(param2')))
%}
end


%% run without inlet, look at 2d coupling with added sinusoid

load('barrier_initial_condition.mat')

%plot([out.x_b_save(:,1:4:40) out.x_s_save(:,1:4:40)])
hold on
cmap = parula(60);
for ii=1:10:60,
patch('XData',[1:500 500:-1:1],'YData',[out.x_s_save(:,ii)' out.x_b_save(500:-1:1,ii)'],'FaceColor',cmap(ii,:),'EdgeColor','black','FaceAlpha',0.5);
end

%% parameter run

name_str = {'storm','slr','wave_height','wave_asym','s_background','w_b_crit','h_b_crit','Qow_max','bb_depth','grain_size','Jmin','a_0','omega','marsh_cover'};

for ii=1,
    
    load(['phase_space3' name_str{ii}])
    figure(ii), hold on
    for kk=1:3,
    for jj=1:length(param1),
        
    out = output{jj,kk};
    if out.x_b_save(50,end)==0,
        Qoverwash(jj) = nan;
        Qinlet(jj) = nan;
    else,
        Qoverwash(jj) = mean(out.Qoverwash);
        Qinlet(jj) = mean(out.Qinlet);
    end
    end
    F{ii}(kk,:) = Qinlet./(Qoverwash+Qinlet)
    end
    plot(param1,F{ii},'-o')
    set(gca,'ylim',[0 1])
    title(name_str{ii})
    clear Qinlet Qoverwash
    
end

a = tight_subplot(4,4,[0.1 0.02]);
for ii=1:1
    load(['phase_space3' name_str{ii}])
    errorbar(a(ii),param1,mean(F{ii}),std(F{ii}))
    xlabel(a(ii),name_str{ii},'interpreter','none')
    set(a(ii),'ylim',[0 0.5])
end
    
%% 

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
