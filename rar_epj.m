% RC_EPJ computes the rotation curves using for a spherical galaxy, a
% Kuzmin disk and an LSB exponential disk using analytic approximations for
% the pattern halo. This code generates figures 7 and 11 in the paper
% "Pattern dark matter and galaxy scaling relations".

% This file is part of PatternDM, a collection of programs for computing
% pattern dark matter halos and associated galactic rotation curves.
% 
% Copyright (C) 2021 by Shankar C. Venkataramani
% <shankar@math.arizona.edu>
% 
% This program is free software; you can redistribute
% it and/or modify it under the terms of the GNU
% General Public License as published by the Free
% Software Foundation; either version 2 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it
% will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General
% Public License along with this program; if not, write
% to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA 02110-1301 USA.

%% Setup and computing the theoretical curves for the RAR plot.

xx = 10.^(-3:0.1:2);

zz = 10.^(-2:0.1:2);

mass = (2*zz-zz.^2.*log((1+zz.^2)./((1+zz).^2))-log((1+zz.^2).*((1+zz).^2)))/4;
% The mass function for the l=0 mode for the Kuzmin disk

g1 = xx./(1-exp(-1*sqrt(xx)));

figure(15);

clf;

% Create axes
axes1 = axes;
hold(axes1,'on');

loglog(xx,xx,'-.','DisplayName','g_{obs}=g_{bar}',...
    'LineWidth',1,'LineStyle','-.',...
    'HandleVisibility','off');

%loglog(xx,sqrt(xx));

%l=sqrt(32*pi);
% l=2;
% 
% loglog(xx,xx+sqrt(xx)-xx.*atan(l./sqrt(xx))/l,...
%     'DisplayName','Compact spherical source','MarkerSize',8,...
%     'Marker','diamond',...
%     'LineStyle','none',...
%     'Color',[1 0 0]);

gam=4;
gbar = 16*(1-1./(1+zz).^(gam-2))./(gam*(gam-3).*zz.^2);
ghal = (2-sqrt(2.)*atan(sqrt(2.)*zz)./zz)./zz;

loglog(gbar,gbar+ghal,...
    'DisplayName','Spherical galaxy','MarkerSize',8,...
    'Marker','diamond',...
    'LineStyle','none',...
    'Color',[1 0 0]);


% loglog(zz./((1+zz.^2).^(3/2)),zz./((1+zz.^2).^(3/2))+zz./(1+zz.^2),...
%     'DisplayName','Critical Kuzmin disk','MarkerSize',8,...
%     'Marker','o',...
%     'LineStyle','-',...
%     'Color',[0 0 1]);

% Uncomment the lines above to get the curve for the "simple" formulae.
loglog(zz./((1+zz.^2).^(3/2)),zz./((1+zz.^2).^(3/2))+mass./(zz.^2),...
    'DisplayName','Critical Kuzmin disk','MarkerSize',8,...
    'Marker','o',...
    'LineStyle','-',...
    'Color',[0 0 1]);

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.001 10]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.001 10]);

box(axes1,'on');
axis(axes1,'square');
% Set the remaining axes properties
set(axes1,'FontSize',16,'XMinorTick','on','XScale','log','YMinorTick','on',...
    'YScale','log');

% Create ylabel
ylabel('$g_{obs}/a_0$',...
    'FontWeight','bold',...
    'FontSize',20,...
    'Interpreter','latex');

% Create xlabel
xlabel('$g_{bar}/a_0$','HorizontalAlignment','center','FontSize',20,...
    'Interpreter','latex');


% RAR fitting function
loglog(xx,g1,'DisplayName','Fitting function Eq.~(50)',...
    'LineWidth',3,...
    'Color',[0 0 0]);




%% LSB exponential disk

% set scales and parameters

a0 = 3600; % in (km/s)^2/kpc

r0 = 1; % r0 in kpc.
v0 = sqrt(a0*r0); %v0 = sqrt(2*pi*G*\Sigmastar*r0)

rmax=7;
A = 1/2;

% exponential disk: Analytic expressions 

vinf = A^(3/4)*v0;
rr = linspace(0,rmax); 

vdisksq = v0^2*A^3*0.5*rr.^2.*(besseli(0,rr/2).*besselk(0,rr/2)-besseli(1,rr/2).*besselk(1,rr/2));
vhalosq = vinf^2*((6+rr)./(3+rr)-6*log(1+rr/3)./rr);
vkepsq = a0*A^3*r0^2./rr; 

adisk = vdisksq./rr;
vrarsq = vdisksq./(1-exp(-sqrt(adisk/a0)));

gdisk = A^3*0.5*zz.*(besseli(0,zz/2).*besselk(0,zz/2)-besseli(1,zz/2).*besselk(1,zz/2));
ghalo = A^(3/2)*((6+zz)./(3+zz)-6*log(1+zz/3)./zz)./zz;


%% Draw figures

figure(15)

loglog(gdisk,gdisk+ghalo,...
    'DisplayName','Exponential disk $(A= 1/2)$','MarkerSize',8,...
    'Marker','*',...
    'Color',[0 0.6 0],...
    'LineStyle','none');


% Create legend
legend1 = legend(axes1,'show');
%set(legend1,'Location','northeastoutside');
set(legend1,'Location','southeast',...
    'Interpreter','latex',...
    'FontSize',14);


% set(legend1,...
%     'Position',[0.447828805852786 0.111184939091917 0.375785725827899 0.322619047619048]);


print -dpng rar_epj.png

figure(13)
clf;

axes2 = axes;
hold(axes2,'on');

plot(rr,sqrt(vdisksq),...
    'DisplayName','$v_{disk}$','Marker','o');
plot(rr,sqrt(vhalosq),...
    'DisplayName','$v_{halo}$','Marker','d','Color',[1 0 0]);
plot(rr,sqrt(vdisksq+vhalosq),...
    'DisplayName','$v_{total}$','Marker','+',...
    'Color',[0.600000023841858 0.200000002980232 0]);

% Create ylabel
ylabel('$v$ in km/s','LineWidth',2,'FontSize',24,...
    'Interpreter','latex');

% Create xlabel
xlabel({'$r$ in kpc'},'FontSize',24,'Interpreter','latex');

box(axes2,'on');
set(axes2,'FontSize',16);

% Create legend
legend2 = legend(axes2,'show');
set(legend2,...
    'Position',[0.45845028252287 0.12089845216069 0.230357279096331 0.28976190884908],...
    'Interpreter','latex',...
    'FontSize',16);

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes2,[0 rmax]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes2,[0 0.6*vinf]);

print -dpng rc_epj.png


