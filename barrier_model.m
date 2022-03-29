function [b_out] = barrier_model(b_struct)


v2struct(b_struct); %retrieve all model parameters back from structure

%dependent variables
u_e_star = u_e./sqrt(g*a0); %equilibrium inlet velocity (non-dimensional)
Vd_max = w_b_crit.*h_b_crit; %max deficit volume m3/m
w_s = R*g*grain_size^2/((18*1e-6)+sqrt(0.75*R*g*(grain_size^3))); %ms-1 church ferguson 2004
phi = 16*e_s*c_s/(15*pi*R*g); %phi from aleja/ashton and trueba/ashton
z0 = 2*wave_height/0.78; % minimum depth of integration (very simple approximation of breaking wave depth based on offshore wave height)
d_sf= 8.9*wave_height; %0.018*wave_height*wave_period*sqrt(g./(R*grain_size)); %depth shoreface m %Hallermeier (1983) or  houston (1995)
k_sf = (3600*24*365)./(d_sf-z0).*(g^(15/4)*wave_height^5*phi*wave_period^(5/2)/(1024*pi^(5/2)*w_s^2) *(4/11*(1/z0^(11/4)-1/(d_sf^(11/4)))));
s_sf_eq = 3*w_s/4/sqrt(d_sf*g)*(5+3*wave_period^2*g/4/(pi^2)/d_sf); %equilibrium shoreface slope
wave_cdf = cumsum(4*[wave_asym*wave_high*ones(wave_climl/4,1);wave_asym*(1-wave_high)*ones(wave_climl/4,1);...
    (1-wave_asym)*(1-wave_high)*ones(wave_climl/4,1);(1-wave_asym)*wave_high*ones(wave_climl/4,1)]./wave_climl);
wave_pdf = 4*[wave_asym*wave_high*ones(wave_climl/4,1);wave_asym*(1-wave_high)*ones(wave_climl/4,1);...
    (1-wave_asym)*(1-wave_high)*ones(wave_climl/4,1);(1-wave_asym)*wave_high*ones(wave_climl/4,1)]./wave_climl;
coast_qs = wave_height.^2.4.*(wave_period.^0.2)*3600*365*24*k .*  (cos(AngArray).^1.2) .* sin(AngArray); %m3/yr
coast_diff = conv(wave_pdf,-(k./(h_b_crit+d_sf)*wave_height^2.4*wave_period^0.2)*365*24*3600 .*  (cos(AngArray).^0.2) .* (1.2*sin(AngArray).^2 - cos(AngArray).^2),'same'); %m2/yr

%timestepping implicit diffusion equation
di = [ny,2:ny,1:ny,1:ny-1,1];
dj = [1,1:ny-1,1:ny,2:ny,ny];

%initial conditions
x_t = (z-d_sf)./s_background + zeros(ny,1); %position shoreface toe m
x_s = rand(ny,1)+d_sf/s_sf_eq+x_t; %position shoreline m
x_b = d_sf/s_sf_eq+w_b_crit+x_t; %position back barrier m
h_b = 2+zeros(ny,1); %height barrier m
barrier_volume = [];
inlet_idx_close_mat = [];
inlet_idx = {};
inlet_idx_mat = [];
inlet_y = zeros(ny,1);

y = 0:dy:(dy*(ny-1)); %alongshore array
t = dt:dt:(dt*nt); %time array

%allocation of output parameters
inlet_nr = zeros(1,length(1:dtsave:nt),'uint16');
inlet_migr = zeros(1,length(1:dtsave:nt),'int16');
inlet_ai = zeros(1,length(1:dtsave:nt),'int32');
inlet_alpha = zeros(1,length(1:dtsave:nt),'single');
inlet_delta = zeros(1,length(1:dtsave:nt),'single');
inlet_beta = zeros(1,length(1:dtsave:nt),'single');
inlet_Qs_in = zeros(1,length(1:dtsave:nt),'single');
inlet_age = cell(nt,1);
Qoverwash = zeros(nt,1,'single');
Qinlet = zeros(nt,1,'single');
c_idx = zeros(ny,1000,'uint8');
bar_strat_x= x_b(1)+1000; %cross-shore location where to record stratigraphy. I guess would be better to do it at one instant in time rather than space?

x_t_save = zeros(ny,length(1:dtsave:nt),'int32'); x_t_save(:,1) = x_t;
x_s_save = zeros(ny,length(1:dtsave:nt),'int32'); x_s_save(:,1) = x_s;
x_b_save = zeros(ny,length(1:dtsave:nt),'int32'); x_b_save(:,1) = x_b;
h_b_save = zeros(ny,length(1:dtsave:nt),'single'); h_b_save(:,1) = h_b;
s_sf_save = zeros(ny,length(1:dtsave:nt),'single'); s_sf_save(:,1) = s_sf_eq;

if plot_on,
    fig = figure('color','white','Name',name,'visible','on',...
        'InvertHardCopy', 'off');
    ax = gca;
    set(ax,'fontsize',15,'Layer','top','XLim',[0 ny*dy/1000],'box','on',...
        'PlotBoxAspectRatio',[2 1 1],'DataAspectRatio',[1 200 1],'nextplot','replacechildren')
    xlabel('Alongshore (km)','fontsize',15)
    ylabel('Cross-shore (m)','fontsize',15)
    if make_gif, set(fig,'Visible','off'), end
end

for i=2:length(t), % years
    
    %show progress for long runs
    if ~mod(i,10000), 
        disp(num2str(i)), end
    
    %sea level
    z = z+(dt*slr); %height of sea level
        
    w = x_b-x_s; %barrier width
    
    d_b = min(bb_depth*ones(size(x_b)),z-(s_background.*x_b)); %basin depth
    s_sf = d_sf./(x_s-x_t); %shoreface slope
    
    
    
    %if the barrier drows, break
    if sum(w<-10)>(ny/2) || any(w<-1000)
        disp('Barrier Drowned')
        break
    end
    
    if barrier_model_on
        %volume deficit
        Vd_b = max(0,(w_b_crit-w).*(h_b+d_b));
        Vd_h = max(0,(h_b_crit-h_b).*w);
        Vd = Vd_b+Vd_h;
        
        %overwash fluxes
        Qow_b = dt.*Qow_max.*Vd_b./max(Vd,Vd_max);
        Qow_h = dt.*Qow_max.*Vd_h./max(Vd,Vd_max);
        Qow = Qow_b+Qow_h;

        %shoreface flux
        Qsf = dt.*k_sf.*(s_sf_eq-s_sf);
        
        %changes
        ff = ((z-s_background.*x_b-d_b)./(z-s_background.*x_b+h_b));
        x_t_dt = (4*Qsf.*(h_b+d_sf)./(d_sf.*(2*h_b+d_sf)))+(2.*dt.*slr./s_sf);
        x_s_dt = 2*Qow./((2*h_b)+d_sf)./(1-ff)-(4*Qsf.*(h_b+d_sf)./(((2*h_b)+d_sf).^2));
        x_b_dt = Qow_b./(h_b+d_b);
        h_b_dt = (Qow_h./w)-(dt.*slr);
        
        %how much q overwash w in total m3/yr
        Qoverwash(i) = sum(dy.*Qow_b./dt);
    else,
        x_t_dt = 0;
        x_s_dt = 0;
        x_b_dt = 0;
        h_b_dt = 0;
    end
    
    if ast_model_on %only alongshore transport calculation to estimate flux into inlets
        
        %simple conv approach
        theta = 180*(atan2([x_s(2:end); x_s(1)]-x_s,dy))/pi;
        %wave direction
        wave_ang = find(wave_cdf>rand,1);
        
        %sed transport this timestep
        Qs = dt.*coast_qs(min(wave_climl,max(1,round(wave_climl-wave_ang-(wave_climl./180.*theta)+1))));
       
    end
    
    
    if inlet_model_on
        %array for changes to back barrier due to flood tidal deltas
        x_b_fld_dt = zeros(ny,1);

        
        %barrier volume is barrier width times height + estimated inlet depth
        barrier_volume = w.*(h_b+2).*sign(min(w,h_b));
        barrier_volume([inlet_idx{:}]) = inf;
        inlet_idx = [inlet_idx num2cell(find(barrier_volume<0)')]; %add drowned barrier to list of inlets
        
        %storm for new inlet every 10 year
        if mod(t(i),10)<(dt/2) && numel(inlet_idx)<inlet_max,
            
            %potential basin length
            if isempty([inlet_idx{:}]),
                basin_length = Jmin+zeros(ny,1);
            else,
                %change%
               basin_length = min(min(Jmin,2* dy*abs(bsxfun(@minus,1:ny,reshape(bsxfun(@plus,[-ny 0 ny],[inlet_idx{:}]'),[],1)))))';
            end
            
            %basin width is simpeler
            basin_width = max(0,z./s_background - x_b); %cross-barrier basin width (m)

            %find new inlets only if its far enough away from existing inlets
            idx = find(basin_length>(Jmin-1));
            [~,new_inlet] = min(barrier_volume(idx));
            new_inlet = idx(new_inlet);
            
            inlet_idx = [inlet_idx {new_inlet}];%add new breach to list of inlets

        end
        
        
        %get rid of duplicates and neighbours
        if ~isempty(inlet_idx)
            
            inlet_idx_mat = [inlet_idx{:}];
            [inlet_all_idx,inlet_all_idx_idx] = sort(inlet_idx_mat);
            %don't try to understand this line.
            inlet_idx_mat(inlet_all_idx_idx(diff([inlet_all_idx(end)-ny inlet_all_idx])<=1)) = [];
            inlet_idx = num2cell(inlet_idx_mat);
        end
        
        %do "fluid mechanics" of inlets
        if ~isempty(inlet_idx),
            
            %sort inlets and find respective tidal prisms
            [inlet_all_idx,inlet_all_idx_idx] = sort(inlet_idx_mat);
            inlet_dist = diff([inlet_all_idx(end)-ny inlet_all_idx inlet_all_idx(1)+ny]);
            
            basin_length = min(Jmin,dy * 0.5 * (inlet_dist(1:end-1) + inlet_dist(2:end))');
            
            %see swart zimmerman
            ah_star = omega0*w(inlet_idx_mat)./sqrt(g*a0);
            c_d = g.*man_n^2./(d_b(inlet_idx_mat).^(1/3));
            
            gam = max(1e-3,inlet_asp.*((omega0.^2).*(1-marsh_cover).^2.*(basin_length(inlet_all_idx_idx).^2).*(basin_width(inlet_idx_mat).^2).*a0./g).^(1/4)./((8/3/pi).*c_d.*w(inlet_idx_mat)));
            a_star_eq = a_star_eq_fun(ah_star,gam,u_e_star);
            u_eq = real(u(a_star_eq,gam,ah_star,a0));
            ai_eq = (omega0.*(1-marsh_cover).*basin_length(inlet_all_idx_idx).*basin_width(inlet_idx_mat).*sqrt(a0./g)).*a_star_eq;
   
            %keep inlet open if velocity is at equilibrium (escoffier). Add
            %margin of 0.05m/s for rounding errors etc
            inlet_close = (u_eq < (u_e-0.05) | isnan(u_eq)) & w(inlet_idx_mat)>0;
            inlet_idx_close_mat = inlet_idx_mat(inlet_close); %close inlets if necessary
            
            %fill inlets back up with alongshore sediment (instantaneously)
            %(in progress)
            %x_s_dt(inlet_idx_close_mat) = x_s_dt(inlet_idx_close_mat)+max(0,w(inlet_idx_close_mat)).*sqrt(ai_eq(inlet_close))*inlet_asp./(d_sf);

            
            %we don't have to think about this one every again!
            inlet_idx(inlet_close) = [];
            inlet_idx_mat(inlet_close) = [];
            ai_eq(inlet_close) = [];
            
            wi_eq = sqrt(ai_eq)./inlet_asp; %calculate width and depths
            di_eq = ai_eq./wi_eq;
            wi_cell = ceil(wi_eq./dy); %get cell widths per inlet
            
        end
        
        
        migr_up = zeros(size(inlet_idx)); %array for inlet migration
        
        for j=1:length(inlet_idx), %inlet morphodynamics per inlet
            
            
            %breach sediment is added to the flood-tidal delta;
            if inlet_idx{j}==new_inlet,
                new_inlet_idx = mod(new_inlet+(1:wi_cell(j))-2,ny)+1;
                x_b_fld_dt(new_inlet_idx) = x_b_fld_dt(new_inlet_idx) + ((h_b(new_inlet)+di_eq(j)).*w(new_inlet))./(d_b(new_inlet));
                
                Qinlet(i) = Qinlet(i) + ((h_b(new_inlet)+d_b(new_inlet)).*w(new_inlet).*wi_cell(j).*dy);
            end
            
            %alongshore flux brought into inlet
            Qs_in(j) = Qs(inlet_idx{j}(1));
            
            %find cells of inlet, updrift barrier, and downdrift barrier
            inlet_idx(j) = {mod(inlet_idx{j}(1)+(1:wi_cell(j))-2,ny)+1};
            inlet_nex(j) = {mod(inlet_idx{j}(end),ny)+1};
            inlet_prv(j) = {mod(inlet_idx{j}(1)-2,ny)+1};

            %find momentum balance of inlet to determine sediment
            %distribution fractions
            Mt = rho_w*u_e.*u_e.*ai_eq(j);
            Mw = rho_w/16*g*wave_height^2.*wi_eq(j);
            I = Mt./Mw.*wi_eq(j)./w(inlet_idx{j}(1)); %test.......
            h_b(inlet_idx{j}) = 0;
            % constrain to not widen
            Ab_prv = w(inlet_prv{j}).*(h_b(inlet_idx{j}(1))+di_eq(j));
            Ab_nex = w(inlet_nex{j}).*(h_b(inlet_nex{j})+di_eq(j));
       
            %do fld delta eq volume
            Vfld = (x_b(inlet_idx{j}(1)) - x_s(inlet_idx{j}(1))+w_b_crit)*wi_eq(j)*d_b(inlet_idx{j}(1));
            Vfld_max = (1e4*(u_e*ai_eq(j)/2/omega0)^0.37); %1e5; %
            
            %add fix to limit unrealistic flood-tidal delta size (based on
            %johnson flood-tidal delta of florida 2006
            if Vfld > Vfld_max
                I = 0.1;
            end
            
            %calculate fractions based on I
            %delta(j) = inlet_fraction(0.05,0.95,3,-3,I);
            %beta(j) = inlet_fraction(0,0.9,10,3,I);
            %beta_r(j) = inlet_fraction(0,0.9,0.9,-3,I);   
            
            delta(j) = inlet_fraction(0,1,3,-3,I);
            beta(j) = inlet_fraction(0,1,10,3,I);
            beta_r(j) = inlet_fraction(0,0.9,0.9,-3,I); 
            
            %{ 
            humans affect inlets?
            delta(j) = 0;
            beta(j) = 1;
            beta_r(j) = 0;
              %}          
            
            alpha(j) = 1-beta(j)-delta(j);
            
            alpha_r(j) = alpha(j)*0.6;
                        
            if Vfld > Vfld_max,
                delta_r(j) = 0;

            else,
                delta_r(j) = ((Ab_nex*(alpha(j)))-(Ab_prv*(beta_r(j))))/Ab_prv;
            end
            
            %use fractions to physically move inlets and fld-tidal detlas.
            fld_delta = abs(Qs_in(j)).*(delta(j)+delta_r(j));%update fld delta, deposit sediment at 0 water depth
            inlet_sink = abs(Qs_in(j)).*(1-beta(j)-beta_r(j));%remove sediment from the shoreface
            
            temp_idx = [inlet_prv{j} inlet_idx{j} inlet_nex{j}];%spread fld tidal delta along one more cell alongshore in both directions
            
            x_b_fld_dt(temp_idx) = x_b_fld_dt(temp_idx)...
                + fld_delta./((numel(temp_idx))*dy)./(h_b(temp_idx)+d_b(temp_idx));
            
            %migrate inlet indices
            migr_up(j) = Qs_in(j).*(alpha_r(j)+alpha(j))./Ab_prv; %migration in m/dt
            migr_dw = Qs_in(j).*(alpha_r(j)+beta_r(j)+delta_r(j))./Ab_nex; %migration in m/dt
            
            %calculate where in the grid cell the inlet is, and add the fractional migration to it
            inlet_y(inlet_idx{j}(1)) = inlet_y(inlet_idx{j}(1))+migr_up(j)./dy;
            
            %how far are the inlets in their gridcell? (or is inlet_y>1 or <0 and should the inlet hop one grid cell?
            migr_int = floor(inlet_y(inlet_idx{j}(1)));
            migr_res = mod(inlet_y(inlet_idx{j}(1)),1);
            
            %reset old grid cell
            inlet_y(inlet_idx{j}(1)) = 0;
            
            %move inlet in gridcell
            inlet_idx{j} = mod(inlet_idx{j} + migr_int - 1,ny) + 1;
            
            inlet_y(inlet_idx{j}(1)) = migr_res;
            
            %how much q flood tidal delta in total
            Qinlet(i) = Qinlet(i)+inlet_sink; %m3 per time step
            
            %add inlet sink to shoreline change
            %x_s_dt(inlet_nex{j}) = x_s_dt(inlet_nex{j})+inlet_sink./(h_b(inlet_nex{j})+d_sf)./dy;

            x_s_dt(temp_idx) = x_s_dt(temp_idx)+inlet_sink./(h_b(temp_idx)+d_sf)./length(temp_idx)./dy;
     
        end
        
        new_inlet = [];
        
        %inlet statistics
        inlet_age{i} = [i*ones(length(inlet_idx_mat),1,'int32'), int32(inlet_idx_mat')]; %fancy lightweight way to keep track of where inlets are in the model
        if mod(i,dtsave)==1,
        inlet_nr(1+fix(i/dtsave)) = length(inlet_idx);
        inlet_migr(1+fix(i/dtsave)) = mean(migr_up./dt);
        if ~isempty(inlet_idx),
            inlet_Qs_in(1+fix(i/dtsave)) = mean(Qs_in);
            inlet_alpha(1+fix(i/dtsave)) = mean(alpha);
            inlet_beta(1+fix(i/dtsave)) = mean(beta);
            inlet_delta(1+fix(i/dtsave)) = mean(delta);
            inlet_ai(1+fix(i/dtsave)) = mean(ai_eq);
            
        end
        end
        
        
    else,
        Qs_in(i)=0;
        delta=0;
        delta_r=0;
        inlet_sink = 0;
        x_b_fld_dt = 0;
    end
    
    %do implicit thing
    if ast_model_on,
        r_ipl = max(0,coast_diff(max(1,min(wave_climl,round(90-theta)))).*dt/2/dy^2); %
        
        dv = [-r_ipl(end); -r_ipl(2:end); 1+2*r_ipl; -r_ipl(1:end-1); -r_ipl(1)];
        A = sparse(di,dj,dv);

        RHS = x_s + r_ipl.*(x_s([2:ny,1]) - 2*x_s + x_s([ny,1:ny-1])) + x_s_dt;
        
        x_s = A\RHS;
    else,
        x_s = x_s + x_s_dt;
    end


    %how are the other moving boundaries changing?
    x_t = x_t + x_t_dt;
    x_b = x_b + x_b_dt + x_b_fld_dt;
    h_b = h_b + h_b_dt;
    
    %to keep track of stratigraphy!
    if sedstrat_on && max(x_b) > bar_strat_x && min(x_s) < bar_strat_x %barrier recording ON
        
        idx_bar = find((x_b-x_b_dt-x_b_fld_dt)<bar_strat_x & (x_b-x_b_fld_dt)>bar_strat_x);
        idx_fld = find((x_b-x_b_dt-x_b_fld_dt)<bar_strat_x & (x_b-x_b_dt)>bar_strat_x);
        idx_h = (x_b-x_b_dt-x_b_fld_dt)>bar_strat_x & x_s<bar_strat_x;
        
        if isempty(inlet_idx),
            idx_inl = [];
        else,
            idx_inl = zeros(size(1,ny));
            for j=1:length(inlet_idx),
                idx_inl = idx_inl | (ismembc(1:ny,inlet_idx{j}) & ((x_s(inlet_nex{j})+w_b_crit)>bar_strat_x)  & ((x_s(inlet_nex{j}))<bar_strat_x));
            end
        end
        
        d_b_strat = round(mean(d_b)/dz_strat);
        z_strat = round(z/dz_strat);
        a0_strat = round(a0/dz_strat);
        if ~isempty(idx_bar),
            for ii=1:length(idx_bar)
                c_idx(idx_bar(ii),d_b_strat:(a0_strat+z_strat)) = 1;
            end
        end
        if ~isempty(idx_fld),
            for ii=1:length(idx_fld)
                c_idx(idx_fld(ii),d_b_strat:(a0_strat+z_strat)) = 2;
            end
        end
        if any(idx_h),
            c_idx(idx_h,z_strat+(a0_strat:(h_b(idx_h)/dz_strat))) = 3;
        end
        if any(idx_inl),
            c_idx(idx_inl,max(1,z_strat+(-round(mean(di_eq/dz_strat)):a0_strat))) = 4;
            c_idx(idx_inl,(a0_strat+z_strat):end) = 0;
        end
        
        if any(inlet_idx_close_mat),
           c_idx(inlet_idx_close_mat,max(1,z_strat+(-round(mean(di_eq/dz_strat)):a0_strat))) = 5;
            c_idx(inlet_idx_close_mat,(a0_strat+z_strat):end) = 0;
        end 
        
    end

    %save variables
    if mod(i,dtsave)==1,
        x_t_save(:,1+fix(i/dtsave)) = x_t;
        x_s_save(:,1+fix(i/dtsave)) = x_s;
        x_b_save(:,1+fix(i/dtsave)) = x_b;
        h_b_save(:,1+fix(i/dtsave)) = h_b;
        s_sf_save(:,1+fix(i/dtsave)) = s_sf;
        
    end
    
    if plot_on && mod(i,100)==1,
        
        if make_gif,
            
            title(['Time: ' num2str(t(i),'%6.0f') 'yr'])
            axis(ax,'auto')
            hx = plot(y./1000,x_s);
            
            
            ymin = min(get(ax,'YLim'));
            
            set(ax,'XLim',get(ax,'XLim'),'YLim',get(ax,'YLim'),'Color',[0 0.5 0]);
            delete(hx);
            
            hh(1) = patch(ax,'XData',[0 0 100 100],'YData',[ymin; z./s_background; z./s_background; ymin]);
            hh(2) = patch(ax,'XData',[0 0 y./1000 100 100],'YData',[ymin; x_b(1); x_b; x_b(end); ymin]);
            hh(3) = patch(ax,'XData',[0 0 y./1000 100 100],'YData',[ymin; x_s(1); x_s; x_s(end); ymin]);
            hh(4) = patch(ax,'XData',[0 0 y./1000 100 100],'YData',[ymin; x_t(1); x_t; x_t(end); ymin]);
            
            
            for j=1:length(inlet_idx),
                if abs(inlet_idx{j}(end)-inlet_idx{j}(1))>(0.5*dy)
                    continue,
                end
                patch([y(inlet_idx{j}(1))./1000 y(inlet_idx{j}(end))./1000 y(inlet_idx{j}(end))./1000 y(inlet_idx{j}(1))./1000],[x_s(inlet_idx{j}(1)) x_s(inlet_idx{j}(end)) x_b(inlet_idx{j}(end)) x_b(inlet_idx{j}(1))],'r','Parent',ax)
            end
            set(hh(1),'FaceColor',[0.08 0.75 1])
            set(hh(2),'FaceColor',[1 1 0])
            set(hh(3),'FaceColor',[0 0 1])
            set(hh(4),'FaceColor',[0 0 0.5])
            lg = legend('Back-barrier','Barrier','Shoreface','Offshore','location','SouthEast');
            set(lg,'color','w')
            
            [imind,cm] = rgb2ind(frame2im(getframe(fig)),256);
            
            if i == 101;
                imwrite(imind,cm,[name '.gif'],'gif', 'Loopcount',inf,'DelayTime',1/15);
            else
                imwrite(imind,cm,[name '.gif'],'gif','WriteMode','append','DelayTime',1/50);
            end
            
        else,
            plot(ax,y./1000,[x_t, x_s, x_b],'-',y([inlet_idx{:} find(barrier_volume<0)'])./1000, x_s([inlet_idx{:} find(barrier_volume<0)']),'o')
            drawnow
        end
        
    end
    
end


b_out = v2struct({'Qoverwash','x_t_save','x_s_save','x_b_save','h_b_save','s_sf_save','z0','d_sf','k_sf','s_sf_eq','fieldNames'});

if sedstrat_on,
    
    %correct for last fill
    [fill_x,fill_y] = find(diff(double([zeros(ny,1), c_idx]),1,2) == -4);
    
    if isempty(di_eq), di_eq = 2; end
    
    for kk = 1:length(fill_x)
        c_idx(fill_x(kk),max(1,(fill_y(kk)-round(mean(di_eq/dz_strat)))):size(c_idx,2)) = 0;
    end
    
    b_out.c_idx = c_idx; 
    b_out.bar_strat_z = z_strat;
    b_out.bar_strat_x = bar_strat_x;
end
if inlet_model_on,
    b_out.inlet_age = cell2mat(inlet_age(~cellfun(@isempty,inlet_age)));
    b_out.inlet_nr = inlet_nr;
    b_out.inlet_migr = inlet_migr;
    b_out.inlet_ai = inlet_ai;
    b_out.inlet_alpha = inlet_alpha;
    b_out.inlet_beta = inlet_beta;
    b_out.inlet_delta = inlet_delta; 
    b_out.inlet_Qs_in = inlet_Qs_in;
    b_out.Qinlet = Qinlet/dt; %put into m3/yr
end

end

