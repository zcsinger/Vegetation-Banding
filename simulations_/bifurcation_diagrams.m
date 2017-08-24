% Generate bifurcation diagram for original model (and possibly more)
% for simple use.

% Zachary Singer, Univerisy of Minnesota - Twin Cities, 8/18/17

%%%%%%%%%%%%%%%%%%%% Specifications and Guide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots original system kinetics with stable and unstable sections.

% 0 = wb^2 - b;

% You can turn on/off the lines.

% Arrows need to be added manually.

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

%%% AESTHETIC PLOTTING CONSTANTS %%%

eq_lines = 'on'; % turn on if you want conserved lines

fontsize = 48;
tick_fontsize = 24;
t_font = 0.86; % what percent of fontsize should title font be (smaller)
linewidth = 4;
stab_color = 'b'; % blue
unstab_color = 'r'; % red

%%% PLOTTING VECTORS AND CONSTANTS %%%

d_b = 0.05;

b_stab_start = 1;
b_stab_stop = 4; % this determines axes length, etc.

b_unstab_start = 1 / b_stab_stop; % will fit perfectly in plot
b_unstab_stop = 1;

b_unstab = b_unstab_start : d_b : b_unstab_stop;
b_stab = b_stab_start : d_b : b_stab_stop;

w_stab = 1 ./ b_stab;
w_unstab = 1 ./ b_unstab;

N = 100;
w_axis = zeros(100, 1);
w_vals = linspace(0, b_stab_stop, N); % N data points on w_axis

c_line = @(c) c - w_vals;

c_start = 0;
d_c = 0.25;
c_stop = 2*b_stab_stop;

%%% PLOTTING %%%

plot(w_axis, w_vals, stab_color, 'LineWidth', linewidth);
hold on
axis([0 b_stab_stop 0 b_stab_stop])
xt = get(gca, 'XTick');
set(gca, 'FontSize', tick_fontsize)
plot(b_stab, w_stab, stab_color, 'LineWidth', linewidth);
plot(b_unstab, w_unstab, unstab_color, 'LineWidth', linewidth);
scatter([1], [1], 'k', 'filled', 'LineWidth', linewidth)
title('Bifurcation diagram for original system', 'FontSize', t_font*fontsize);
xlabel('Biomass $$b$$', 'FontSize', fontsize);
ylabel('Water $$w$$', 'FontSize', fontsize);


if eq_lines == 'on'
    
    for c = c_start : d_c : c_stop
        c_val = c_line(c);
        plot(w_vals, c_val, 'k');
    end
    
end

