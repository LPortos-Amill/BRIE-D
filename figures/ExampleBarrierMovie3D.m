%ExampleBarrierMovie
clr
name = 'DataS1_ExampleBarrierPlot';
name = 'drowning_inlet';
load([dropbox filesep 'BarrierModel' filesep name '.mat'])
out = output{2,6};
b_struct.ny = 1000;
b_struct.nt = 1e5;
b_struct.dtsave = 1000;
b_struct.dt = 0.05;
b_struct.dy = 100;

step = 10;

x_s = out.x_s_save./step;
x_t = out.x_t_save./step;
x_b = out.x_b_save./step;

t = (1:b_struct.dtsave:b_struct.nt)*b_struct.dt;
z = b_struct.z+(t*b_struct.slr);

x_grid = 0:step:(max(z)/b_struct.s_background);
y_grid = (1:b_struct.ny)*b_struct.dy;

fig = figure('color','white','visible','off',...
    'InvertHardCopy', 'off','Position',[100 100 1200 800]);
ax = axes('fontsize',9,'Layer','top','box','on','nextplot','replacechildren');

set(ax,'DataAspectRatio',[1 1 10],'View',[-40 30]),
demcmap([-10 5],150); %hold on
l1 = xlabel('\leftarrow Ocean | Land (km) \rightarrow','FontSize',10); %,'position',[12.0639   -2.8840  -17.8682]);
l2 = ylabel('Alongshore (km)','FontSize',10); %,'position',[-1.8547   25.2491  -16.1331]);
l3 = zlabel('Elevation (m)','FontSize',10); %,'position',[-1.9708   51.5122   10.0000]);
%align_axislabel([],ax)
annotation('textbox', [0.05, 0, 0.1, 0.1], 'string', 'BRIE model, Nienhuis & Lorenzo-Trueba, GRL 2020','edgecolor','none','fontsize',10)
cmap = demcmap([-10 5],150);
cmap = [cmap(1:100,:); cmap(round(125:-0.5:101),:)];
colormap(cmap)
cbar = colorbar('Location','east','position',[0.75 0.55 0.02 0.25],'Fontsize',8);

for tt=1:length(t),
    tt
    title(['Time: ' num2str(t(tt),'%6.0f') ' yr          Cumulative sea-level rise: ' num2str(z(tt)-10,'%1.1f') ' m'],'FontSize',11)
    z_grid = repmat(x_grid*b_struct.s_background,[length(y_grid),1]);
    %grid
    for yy = 1:length(y_grid),
        idx = x_t(yy,tt):x_s(yy,tt);
        z_grid(yy,idx) = z_grid(yy,idx(1)) + out.s_sf_save(yy,tt)*(1:length(idx)).*step;
        
        idx = x_s(yy,tt):x_b(yy,tt);
        z_grid(yy,idx) = out.h_b_save(yy,tt)+z(tt);
        
        idx = x_b(yy,tt):round((z(tt)-b_struct.bb_depth)/b_struct.s_background./step);
        z_grid(yy,idx) = z(tt)-b_struct.bb_depth;
        
    end
    idx_grid = [1:step:(max(1,x_t(1,tt)-(1000/step))),...
        (max(1,x_t(1,tt)-(1000/step))):round(z(tt)/b_struct.s_background./step),...
        round(z(tt)/b_struct.s_background./step):step:length(x_grid)];
    
    %inlets
    inl_idx = double(out.inlet_age(out.inlet_age(:,1)==(tt*b_struct.dtsave),2));
    inl_w = sqrt(double(out.inlet_ai(tt)))/b_struct.inlet_asp; w_idx = round(inl_w./b_struct.dy);
    inl_d = sqrt(double(out.inlet_ai(tt)))*b_struct.inlet_asp;
    
    for j=1:length(inl_idx),
        
        
        fld_all = max(1,min(b_struct.ny,inl_idx(j)+round(-(w_idx)+(0:(2*w_idx)))));
        [x_all,y_all] = colon_all(min(x_b(fld_all,tt)).*ones(length(fld_all),1,'int32'),x_b(fld_all,tt),fld_all,step);
        
        z_grid(y_all,x_all) = z(tt)-0.5;
        
        inl_all = max(1,min(b_struct.ny,inl_idx(j)+round(-(w_idx/2))+(0:w_idx))); %
                
        [x_all,y_all] = colon_all(x_s(inl_all,tt)-50,x_b(inl_all,tt),inl_all,step);
        
        z_grid(y_all,x_all) = z(tt)-inl_d;
       
    end
    %figure, imagesc(z_grid-z(tt),[-3 3])
    set(ax,'CLim',[-10 5]+z(tt)); set(cbar,'Ticks',[-10 -5 0 5]+z(tt),'TickLabels',{'-10','-5','0','5'});
    cbar.Label.String = 'Meters above sea level';
    cbar.Label.FontSize = 10;
    
    surf(ax,x_grid(idx_grid)./1000,y_grid./1000,z_grid(:,idx_grid),'EdgeColor','none','FaceColor','interp')
    %drawnow
    %%{
    if tt == 1
        
        [A,map] = rgb2ind(frame2im(getframe(fig)),256);
        imwrite(A,map,[name '3D.gif'],'gif','LoopCount',Inf,'DelayTime',1/10);
    else
        [A,map] = rgb2ind(frame2im(getframe(fig)),256);
        imwrite(A,map,[name '3D.gif'],'gif','WriteMode','append','DelayTime',1/20);
    end
    %}

    
end