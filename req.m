% REQ computes the equivalent radius for spherical galaxies for various
% choices of the exponent gamma. This code generates figure 9 in the paper
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

%% Setup 

gam_vals = [3+1./(5:-1:2) 4:13];
ngam = length(gam_vals);

rhalf = 2.^(1./(gam_vals-3))-1;
re = zeros(1,ngam);

for i=1:ngam
    
    gam = gam_vals(i);
    
    f = chebfun2(@(re,q) (cos(q)./(cos(q)+re+10^(-6))).^(gam-2),[0 rhalf(i) -pi/2 pi/2],'splitting','on');
    
    g=sum(f);
    
    h=cumsum(g);
    
    re(i) = roots(h-1/(gam-3));
    
end

figure(5)
clf;
axes1 = gca;

plot(-1.6:0.05:2.4,log(2.^(exp(1.6:-0.05:-2.4))-1),'linewidth',3);

hold on;

plot(log(gam_vals-3),log(re),'r*','MarkerSize',8);


% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[-2 2.5]);
box(axes1,'on');
% Set the remaining axes properties
%set(axes1,'FontSize',14,'PlotBoxAspectRatio',[1 0.65 0.5],'XAxisLocation',...
%   'origin','YAxisLocation','origin');
set(axes1,'FontSize',14,'XAxisLocation','origin','YAxisLocation','origin');
%'PlotBoxAspectRatio',[1 0.65 0.5],'XAxisLocation','origin','YAxisLocation','origin');

% Create ylabel
ylabel('$\ln(r_e/r_0)$','interpreter','latex','FontSize',18);

% Create xlabel
xlabel('$\ln(\gamma-3)$','interpreter','latex','FontSize',18);

print -dpng req.png