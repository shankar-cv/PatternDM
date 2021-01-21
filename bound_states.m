% bound_states computes the first and second energy level for 1d gaussian
% wells. This code generates figure 12 in the paper
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

pot_vals = 0.05:0.05:1.25;
npot = length(pot_vals);

x = chebfun('x',[-25,25]);
V = -exp(-0.5*x.^2);

energies = zeros(2,npot);

for i=1:npot
    
    W = pot_vals(i)*V;   
    energies(:,i) = quantumstates(W,2,1);
    
end

figure(2)

clf;

plot(pot_vals,-energies);

pot_vals = 1.3:0.05:1.6;
npot = length(pot_vals);

x = chebfun('x',[-100,100]);
V = -exp(-0.5*x.^2);

energies2 = zeros(2,npot);

for i=1:npot
    
    W = pot_vals(i)*V;   
    energies2(:,i) = quantumstates(W,2,1);
    
end

pot_vals = 1.65:0.05:2;
npot = length(pot_vals);

x = chebfun('x',[-25,25]);
V = -exp(-0.5*x.^2);

energies3 = zeros(2,npot);

for i=1:npot
    
    W = pot_vals(i)*V;   
    energies3(:,i) = quantumstates(W,2,1);
    
end

vv = 1.3:0.05:2;
e2 = [-energies2(2,:) -energies3(2,:)];

figure(3);
clf;
axes2=gca;


plot(vv,e2,'linewidth',3,'DisplayName','Energy level');
hold on;
plot(vv,0.154*(vv-1.3).^2,'-.','linewidth',3,'DisplayName','Quadratic fit');

box(axes2,'on');
% Set the remaining axes properties
set(axes2,'FontSize',16,'PlotBoxAspectRatio',[1 0.65 0.5]);

% Create ylabel
ylabel('$\lambda_1=a_0^2 k_0^2$','interpreter','latex','FontSize',20,'FontWeight','bold');

% Create xlabel
xlabel('$\eta = a^2 V_0$','interpreter','latex','FontSize',20,'FontWeight','bold');

% Create legend
legend1 = legend(axes2,'show');
%set(legend1,'Location','northeastoutside');
set(legend1,'Location','southeast',...
    'Interpreter','latex',...
    'FontSize',16);

print -dpng bound-states.png