% Simulation of modified Klausmeier model of two-species diffusion-reaction
% Looking for heteroclinic orbit opposite to that of 
% veg_dir_sim_h1.m

% Zachary Singer, University of Minnesota - Twin Cities, 8/8/17

%%%%%%%%%%%%%%%%%%%% Specifications and Guide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created to verify relationship between water front height w_plus
% and front speed s. 

% Relies on functions in veg_model_f.m and veg_model_j.m

% How to use:

% There are modes you can turn "on" or "off:"
% plotting  - see plots of space x vs. biomass b and space x vs. water w
%
% speed     - find front speed s by using time needed to cross region
%           - You may want to adjust the region where speed is tested 
%             by changing the value t_size from floor(L/4) to floor(L/6) 
%             (accuracy increases value decreases but could be too long)
%
% speed_alt - find front speed s by instead finding the height of the water
%           - interface w_minus, can use explicit formula. Also because
%           - finding speed in same way as veg_dir_sim_h1 doesn't work
%           - (interface is high on left as well as right)
%
% write_data - write results to a text file with the first row being the
%              range of w_plus values desired, second row s values.
%              Change the filename variable to a suitable name.

% Suggested Ranges of Variables:
% c : 1 to 5
% w_minus : 2 to 8
% tmax : 10 to 20

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Original system with traveling wave ansatz: 
% b_t = d_b*b_xx + w*b^2 - mb + sb_x
% w_t =          - w*b^2 + mb + (s+c)w_x

% Scaled out : d_b = 1, m = 1 (for simplicity)

% Analyzed model: 
% b_t = b_xx + w*b^2 - b 
% w_t =      - w*b^2 + b + cw_x

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

% Modes (see above)

plotting = "on";
speed = "on";   % found by measuring time for front to pass section
speed_alt = "off"; % have only one of speed and speed_alt on.
                   % found by finding b_- and w_- heights, uses formula
write_data = "off";
filename = "dirsim_het2_date_01.txt";

%%% AESTHETIC PLOTTING CONSTANTS %%%

fontsize = 36;
tick_fontsize = 20;
linewidth = 4;

pause_time = 7.5; % pause to save figure at certain time
tmax = 50;        % how long to go for

%%% IMPORTANT CONSTANTS TO CHANGE / KEEP TRACK OF %%%

c = 2;  % advective constant
b_plus = 4;


% Range of values

b_plus_s = b_plus;        % starting point for b_plus
d_b = 1;
b_plus_e = b_plus_s; % if you want range of b_plus values
b_range = b_plus_s : d_b : b_plus_e;

b_data = b_range; % b_plus values stored here
s_data = [];      % front speeds stored here

% Space Step Constants

L = 100;             % should be 50 or greater
dx = 0.1;
x = (dx: dx: L)';    % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2);        % how many time steps

% Time Step Constants

dt = .05;
tspan = [0 dt];

% Speed Time Constants

t_eps = 0.001;    % threshold for front when calculating speed s
t_start = 0;    % set to 0 for now
t_stop = 0;     
s = 0;          % placeholder

speed_width = 0.2;  % what percent of domain to measure speed over
speed_start = 0.1;  % what percent of domain to wait before starting

%%%%%% !!! How far between speed measurements !!! %%%%%

t_size = floor(L*speed_width) / dx; 

%%%%%% !!! How far before starting measurements !!! %%%%%

t_thres = floor(L*speed_start) / dx;  

% Derivative Matrices

e = ones(N, 1);

D1 = spdiags([-e e], 0 : 1, N, N) / dx; % first derivative, upwind

D2 = spdiags([e -2*e e], -1 : 1, N, N) / dx^2; % second deriv., central

%%%%%%%%%%%%%%%%%%% Loop over w_plus values: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b_plus = b_range
    
w_plus = 1 / b_plus; % consequence of equilibrium

% Kinetics (found in veg_model_f.m and veg_model_j.m)

F = @(t, V) veg_model_f2(t, V, D1, D2, N, dx, c, b_plus);
J = @(t, V) veg_model_j(t, V, D1, D2, N, c);

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'Jacobian', J);

%%%%% INITIAL INTERFACES FOR BIOMASS AND WATER %%%%%

const = 5; % steepness of interface in middle

b = ((b_plus / 2) * tanh(const*(x - floor(L/4)))) + b_plus / 2 ;
w = (-(b_plus / 2) * (tanh(const*(x - floor(L/4))))) + b_plus / 2 + 1 / b_plus;

U = [b; w]; % how we store biomass and water

if plotting == "on"

%%% INITIAL PLOTS FOR INTERFACES %%%

subplot(2, 1, 1)
plot(x, w, 'b', 'LineWidth', linewidth)  
axis([0 L 0 2*b_plus])
xt = get(gca, 'XTick');
set(gca, 'FontSize', tick_fontsize)
xlabel('Space $$x$$', 'FontSize', fontsize)
ylabel('Water $$w$$', 'FontSize', fontsize)
% title(['Initial Space $$x$$ vs. Water $$w$$ with $$c =$$ ' num2str(c) ...
%         ', $$b_+ =$$ ' num2str(b_plus) ' and time $$t = 0$$'], ...
%         'FontSize', fontsize) % more descriptive title
 title(['Space $$x$$ vs. Water $$w$$ at $$t = 0$$'], 'FontSize', fontsize)

subplot(2, 1, 2)
plot(x, b, 'g', 'LineWidth', linewidth)
axis([0 L 0 2*b_plus]) % Might have to increase y-axis
xt = get(gca, 'XTick');
set(gca, 'FontSize', tick_fontsize)
xlabel('Space $$x$$', 'FontSize', fontsize)
ylabel('Biomass $$b$$', 'FontSize', fontsize)
% title(['Initial Space $$x$$ vs. Biomass $$b$$ with $$c =$$ ' num2str(c) ...
%         ', $$b_+ =$$ ' num2str(b_plus) ' and time $$t = 0$$'], ...
%         'FontSize', fontsize) % more descriptive title
 title(['Space $$x$$ vs. Biomass $$b$$ at $$t = 0$$'], 'FontSize', fontsize)

drawnow
pause

end

%%%%%%%%%%%%%%%%%%%%%% Main Time Stepping:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1 : tmax / dt
      
    [T, UALL] = ode15s(F, tspan + j*dt, U, options);
    U = UALL(end,:)';
   
    b = U(1: N);
    w = U(N+1: 2*N);
    
    %%%%%% SPEED CALCULATIONS (TURN speed = "on" above) %%%%%%
    
    if speed == "on"
        
    if (t_start ~= 0 && t_stop ~= 0)
        continue
    end
        
    % speed calculation (now about b)
    b_sp = find(b > b_plus - t_eps);
    bf_size = size(b_sp);
    bf_len = bf_size(1);
    
    if t_start == 0 && bf_len > 0
        if b_sp(1) < (floor(N/4) + t_thres)
            b_sp(end);
            t_start = j*dt;
        end
    elseif t_stop == 0 && bf_len > 0 
        if b_sp(1) > (floor(N/4) + (t_thres + t_size)) 
            t_stop = j*dt;
        end
    end
    % if we have passed the desired distance find speed s
    if t_start ~= 0 && t_stop ~= 0
        s = (t_size*dx) / (t_stop - t_start);
        s_data = [s_data s];
    end
    
    end
 
    %%%%% PLOTTING (turn plot = "on" above) %%%%%
    
    if plotting == "on"
    
    subplot(2, 1, 1)
    plot(x, w, 'b', 'LineWidth', linewidth)    
    axis([0 L 0 2*b_plus])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', tick_fontsize)
    xlabel('Space $$x$$', 'FontSize', fontsize)
    ylabel('Water $$w$$', 'FontSize', fontsize)
%     title(['Space $$x$$ vs. Water $$w$$ with $$c =$$ ' num2str(c) ...
%             ', $$b_+ =$$ ' num2str(b_plus) ...
%             ', time $$t =$$ ' num2str(j*dt) ', and speed = ' num2str(s)], ...
%             'FontSize', fontsize)
    title(['Space $$x$$ vs. Water $$w$$ at $$t =$$ ' num2str(j*dt)], ...
            'FontSize', fontsize)

    
    subplot(2, 1, 2)
    plot(x, b, 'g', 'LineWidth', linewidth)
    axis([0 L 0 2*b_plus])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', tick_fontsize) 
    xlabel('Space $$x$$', 'FontSize', fontsize)
    ylabel('Biomass $$b$$', 'FontSize', fontsize)
%     title(['Space $$x$$ vs. Biomass $$b$$ with $$c =$$ ' num2str(c) ...
%             ', $$b_+ =$$ ' num2str(b_plus) ', time $$t =$$ ' ...
%              num2str(j*dt) ', and speed = ' num2str(s)], ...
%              'FontSize', fontsize)
    title(['Space $$x$$ vs. Biomass $$b$$ at $$t =$$ ' num2str(j*dt)], ...
            'FontSize', fontsize)
        
    drawnow
    % pause
    
    if j*dt == pause_time
        pause
    end
    
    end
    
end

%%% Alternate way of calculating speed by looking at w_minus height %%%
% Explicit formula: s = c*(w_minus - w_plus) / (b_plus + w_plus - w_minus)

if speed_alt == 'on'
   w = U(N+1 : 2*N);
   w_minus = w(floor(N/5)); % look at flat section
   w_plus = 1 / b_plus;
   
   s = c * (w_minus - w_plus) / (b_plus + w_plus - w_minus);
   s_data = [s_data s];
end

title(['Space $$x$$ vs. Biomass $$b$$ at $$t =$$ ' num2str(tmax)], ...
        'FontSize', fontsize)
subplot(2, 1, 1)
title(['Space $$x$$ vs. Water $$w$$ at $$t =$$ ' num2str(tmax)], ...
        'FontSize', fontsize)


end

%%%%%% WRITE DATA (turn write_data = "on" above) %%%%%

if write_data == "on"

results = [b_data; s_data];
dlmwrite(filename, results);

end
