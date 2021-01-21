% EPJ_FIGS generates figures 3,4 and 8 in the paper
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

%% Figure 8. eta vs zeta

xx = 10.^(-2:0.1:2);
jj = max(1,xx);
yy = (pi^2/4+jj)./(xx.*sqrt(jj));

figure(8)
clf;
axes1=gca;

loglog(xx,yy,'linewidth',3);

% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.1 400]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14,'PlotBoxAspectRatio',[1 0.65 0.5]);

% Create ylabel
ylabel('$\eta$','interpreter','latex','FontSize',20,'FontWeight','bold');

% Create xlabel
xlabel('$\zeta$','interpreter','latex','FontSize',20,'FontWeight','bold');

print -dpng conditions.png

%% Figure 3. Wavevector for target pattern

N = chebop(0.05,10);
N.op = @(x,k) diff(k,2)+2*diff(k)./x - 2*k*(k.^2+1./x.^2 -1);

N.lbc=0.711382*0.05; %magic number from shooting method on mathematica
N.rbc=0.995; % Far field behavior is 1-0.5*xi^(-2)

[k,info] = solvebvp(N,0);

figure(3);
clf;
axes5=gca;


loglog(k,'linewidth',3);

% Uncomment the following line to preserve the limits of the axes
xlim(axes5,[0.05,10]);
ylim(axes5,[0.03 1.3]);
box(axes5,'on');
% Set the remaining axes properties
set(axes5,'FontSize',14,'PlotBoxAspectRatio',[1 0.65 0.5]);

% Create ylabel
ylabel('$k(\xi)/k_0$','interpreter','latex','FontSize',20,'FontWeight','bold');

% Create xlabel
xlabel('$\xi = k_0 R$','interpreter','latex','FontSize',20,'FontWeight','bold');

print -dpng k-halo.png

%% Figure 4. Effective halo density

figure(4);
clf;
axes3=gca;

r = chebfun({@(x) x},[0.05 10]);
rhop = 2*(k.^4 + (diff(k)+2*k./r).^2 - 1);

loglog(rhop,'linewidth',3,'DisplayName','$\rho_P/(\Sigma^* k_0)$');
hold on;
loglog(8./(1+2*r.^2),'--','linewidth',3);

% Uncomment the following line to preserve the limits of the axes
xlim(axes3,[0.05,10]);
ylim(axes3,[0.03 10]);
box(axes3,'on');
% Set the remaining axes properties
set(axes3,'FontSize',14,'PlotBoxAspectRatio',[1 0.65 0.5]);

% Create ylabel
ylabel('$\rho_P/(\Sigma^* k_0)$','interpreter','latex','FontSize',20,'FontWeight','bold');

% Create xlabel
xlabel('$\xi = k_0 R$','interpreter','latex','FontSize',20,'FontWeight','bold');

print -dpng rho-halo.png
