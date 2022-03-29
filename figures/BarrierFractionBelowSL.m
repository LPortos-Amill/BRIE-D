inlet_asp = 100;
w_barrier = 500; %width barrier at inlet (m)
basin_width = 10000; %m
basin_length = 10000; %m
a0 = 0:0.1:1; %amplitude of tide (m)
omega0 = 1.4e-4; %tidal frequency (rad s-1);
gamma = sqrt(0.005); %aspect ratio inlet
g = 9.81;
man_n = 0.05; %sm^-(1/3) manning n (vegetated channel)
u_e = 1; %ms-1 inlet equilibrium velocity (see swart/zimmerman)
u_e_star = u_e./sqrt(g*a0);%new explicit relationship between boundary conditions and inlet area
f_marsh = 0.5;
ah_star = omega0*w_barrier./sqrt(g*a0);
c_d = g.*man_n^2./(2.^(1/3));
gam = gamma.*((omega0.^2).*(1-f_marsh).*(basin_length.^2).*(basin_width.^2).*a0./g).^(1/4)./((8/3/pi).*c_d.*w_barrier);
%
a_star = real(((2*ah_star)./3 + (2^(2/3)*((18*ah_star.*gam.^2 - 27*u_e_star.^4 - 2*ah_star.^3.*gam.^2.*u_e_star.^2 +...
    3*3^(1/2)*gam.^2.*u_e_star.^2.*(-(4*ah_star.^4.*gam.^4.*u_e_star.^4 - 4*ah_star.^3.*gam.^2.*u_e_star.^8 -...
    8*ah_star.^2.*gam.^4.*u_e_star.^2 + 36*ah_star.*gam.^2.*u_e_star.^6 + 4.*gam.^4 - 27.*u_e_star.^10)./...
    (gam.^4.*u_e_star.^6)).^(1/2))./(gam.^2.*u_e_star.^2)).^(1/3))./6 + (2^(1/3).*(ah_star.^2.*u_e_star.^2 + 3))./...
    (3.*u_e_star.^2.*((18*ah_star.*gam.^2 - 27.*u_e_star.^4 - 2*ah_star.^3.*gam.^2.*u_e_star.^2 + 3*3.^(1/2).*gam.^2.*u_e_star.^2.*...
    (-(4*ah_star.^4.*gam.^4.*u_e_star.^4 - 4.*ah_star.^3.*gam.^2.*u_e_star.^8 - 8.*ah_star.^2.*gam.^4.*u_e_star.^2 + 36*ah_star.*gam.^2.*u_e_star.^6 +...
    4.*gam.^4 - 27.*u_e_star.^10)./(gam.^4.*u_e_star.^6)).^(1/2))./(gam.^2.*u_e_star.^2)).^(1/3))));

u = sqrt(g*a0).*sqrt(gam./2.*sqrt(a_star).*((-gam.*sqrt(a_star).*((a_star-ah_star).^2))+sqrt((gam.^2).*a_star.*((a_star-ah_star).^4)+4)));

%inlet cross-sectional area
a = omega0.*basin_width.*basin_length.*(1-f_marsh).*sqrt(a0/g).*a_star;

%fraction of tidal prism compared to potential tidal prism (area*2*a0)
Pfrac = (2.*a.*u./omega0)./(2.*basin_width.*(1-f_marsh).*basin_length.*a0);

%inlet width
w_inlet = sqrt(a)./gamma;


w_inlet_frac = w_inlet./basin_length;

plot(a0,w_inlet_frac.*100)
xlabel('Tidal Amplitude (m)')
ylabel('Barrier fraction below SL (%)')
saveas(gcf,'BarrierFractionBelowSL.png')