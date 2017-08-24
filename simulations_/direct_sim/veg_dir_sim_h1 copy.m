% Simulation of modified Klausmeier model of two-species diffusion-reaction

% Zachary Singer, University of Minnesota - Twin Cities, 8/3/17

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
% write_data - write results to a text file with the first row being the
%              range of w_plus values desired, second row s values.
%              Change the filename variable to a suitable name.
%              NOTE: Does not save biomass or water interface heights...

% Suggested Ranges of Variables:
% c : 1 to 5
% w_plus : 2 to 8
% tmax : 10 to 20

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyzed model: 

% b_t = b_xx + w*b^2 - b 
% w_t = cw_x - w*b^2 + b 

% Key parameters:

% c : advection speed of water
% w_plus : initial height of water interface on right

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up window, interpreter, figure for plotting

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

%%% MODES (see above) %%%

plotting = "on";
speed = "off";
write_data = "off";
filename = "dirsim_het1_date_01.txt";

%%% AESTHETIC PLOTTING CONSTANTS %%%

fontsize = 36;
tick_fontsize = 20;
linewidth = 4;

pause_time = 5; % pause here to save figure

%%%%% IMPORTANT CONSTANTS TO PLAY WITH %%%%%

c = 2;          % advective constant
w_plus_s = 4;   % initial water height
b_minus = 8;    % this might be a relevant height -> 5.83176

tmax = 15; % how long time is

speed_width = 0.2;  % what percent of domain to measure speed over
speed_start = 0.1;  % what percent of domain to wait before starting

%%%%% OTHER CONSTANTS %%%%%

w_minus = 1 / b_minus;

% Right now we just do one w_plus but can make range

d_w = 0.1;
w_plus_e = w_plus_s; 
w_range = w_plus_s : d_w : w_plus_e; 

w_data = w_range; % data about w_plus values used stored here
s_data = [];

% Space Step Constants

L = 100;             % should be 50 or greater
dx = 0.1;
x = (dx: dx: L)';    % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2);        % how many time steps

% Time Step Constants

dt = .05;   % time step
tspan = [0 dt];

% Speed Time Constants

% we set most of these to zero for the time being (part of if statements)
t_eps = 0.1;    % threshold for front when calculating speed s
t_start = 0;    
t_stop = 0;    
s = 0;          

t_size = floor(L*speed_width) / dx;     % how far between speed meas.
t_thres = floor(L*speed_start) / dx;    % how far before starting speed meas.

% Derivative Matrices

e = ones(N, 1);

D1 = spdiags([-e e], 0 : 1, N, N) / dx; % first derivative, upwind

D2 = spdiags([e -2*e e], -1 : 1, N, N) / dx^2; % second deriv., central

%%%%%%%%%%%%%%%%%%% Loop over w_plus values: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for w_plus = w_range

% Kinetics (found in veg_model_f.m and veg_model_j.m)

F = @(t, V) veg_model_f(t, V, D1, D2, N, dx, c, w_plus, b_minus);
J = @(t, V) veg_model_j(t, V, D1, D2, N, c);

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'Jacobian', J);

%%% INTERFACES FOR BIOMASS b AND WATER w %%%

const = 5; % steepness of interface in middle

b = -(b_minus / 2) * (tanh(const*(x - floor(L/4))) + 1) + b_minus ;
w = ((w_plus - w_minus) / 2) * (tanh(const*(x - floor(L/4))) + 1) + w_minus;

U = [b; w];

if plotting == "on"

% Initial Plot to view Interfaces

subplot(2, 1, 1)
plot(x, w, 'b', 'LineWidth', linewidth)    
axis([0 L 0 2.5*w_plus])
xt = get(gca, 'XTick');
set(gca, 'FontSize', tick_fontsize)  % set tick size
xlabel('Space $$x$$', 'FontSize', fontsize)
ylabel('Water $$w$$', 'FontSize', fontsize)
title(['Space $$x$$ vs. Water $$w$$ at $$t = 0$$'], 'FontSize', fontsize)

% More detailed title

% title(['Space $$x$$ vs. Water $$w$$ with $$c =$$ ' num2str(c) ...
%         ', $$w_+ =$$ ' num2str(w_plus) ...
%         ', time $$t = 0$$'], ...
%         'FontSize', fontsize)

subplot(2,1,2)
plot(x, b, 'g', 'LineWidth', linewidth)
axis([0 L 0 2.5*w_plus])
xt = get(gca, 'XTick');
set(gca, 'FontSize', tick_fontsize)
xlabel('Space $$x$$', 'FontSize', fontsize)
ylabel('Biomass $$b$$', 'FontSize', fontsize)
title(['Space $$x$$ vs. Biomass $$b$$ at $$t = 0$$'], 'FontSize', fontsize)

% title(['Space $$x$$ vs. Biomass $$b$$ with $$c =$$ ' num2str(c) ...
%         ', $$w_+ =$$ ' num2str(w_plus) ', time $$t = 0$$'], ...
%             'FontSize', fontsize)

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
    
    % OVERVIEW: Speed is found by looking for first point in w vector
    %           where value is at w_plus value. If this exists, we track
    %           the first time it passes our initial point + threshold,
    %           and look for time it passes our init. + thres. + size
    
    if speed == "on"
        
    if (t_start ~= 0 && t_stop ~= 0) % already found speed, so stop
        break
    end
        
    % Speed calculation
    
    w_sp = find(w > w_plus - t_eps);    % indices of water values above thres
    wf_size = size(w_sp);
    wf_len = wf_size(1);
    
    if t_start == 0 && wf_len > 0       % at least one value and first time
        if w_sp(1) > (floor(N/4) + t_thres)
            t_start = j*dt; % start time
        end
        
    elseif t_stop == 0 && wf_len > 0 
        if w_sp(1) > (floor(N/4) + t_thres + t_size) 
            t_stop = j*dt;  % end time
        end
        
    end
    
    % If we have passed the desired distance, find speed s
    
    if t_start ~= 0 && t_stop ~= 0
        s = (t_size*dx) / (t_stop - t_start);
        s_data = [s_data s];
    end
    
    end
     
    %%%%% PLOTTING (turn plot = "on" above) %%%%%
    
    if plotting == "on"
    
    subplot(2, 1, 1)
    plot(x, w, 'b', 'LineWidth', linewidth)    
    axis([0 L 0 2.5*w_plus])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', tick_fontsize)
    xlabel('Space $$x$$', 'FontSize', fontsize)
    ylabel('Water $$w$$', 'FontSize', fontsize)
    title(['Space $$x$$ vs. Water $$w$$ at $$t =$$ ' num2str(j*dt)], ...
            'FontSize', fontsize)
        
    %     title(['Space $$x$$ vs. Water $$w$$ with $$c =$$ ' num2str(c) ...
    %         ', $$w_+ =$$ ' num2str(w_plus), time $$t =$$ '  ...
    %          num2str(j*dt) ', and speed = ' num2str(s)], ...
    %         'FontSize', fontsize)    
    
    subplot(2, 1, 2)
    plot(x, b, 'g', 'LineWidth', linewidth)
    axis([0 L 0 2.5*w_plus])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', tick_fontsize) 
    xlabel('Space $$x$$', 'FontSize', fontsize)
    ylabel('Biomass $$b$$', 'FontSize', fontsize)
    title(['Space $$x$$ vs. Biomass $$b$$ at $$t =$$ ' num2str(j*dt)], ...
            'FontSize', fontsize)
        
    %     title(['Space $$x$$ vs. Biomass $$b$$ with $$c =$$ ' num2str(c) ...
    %             ', $$w_+ =$$ ' num2str(w_plus) ', time $$t =$$ ' ...
    %              num2str(j*dt) ', and speed = ' num2str(s)], ...
    %              'FontSize', fontsize)    
        
    drawnow
    
    % Pause if we reach pause_time
    
    if j*dt == pause_time
        pause
    end
    
    end
    
end

end

%%%%%% WRITE DATA (turn write_data = "on" above) %%%%%

if write_data == "on"

results = [w_data; s_data];
dlmwrite(filename, results);

end

if speed == 'off'

% Set title to end time
    
title(['Space $$x$$ vs. Biomass $$b$$ at $$t =$$ ' num2str(j*dt)], ...
        'FontSize', fontsize)
subplot(2, 1, 1)
title(['Space $$x$$ vs. Water $$w$$ at $$t =$$ ' num2str(j*dt)], ...
        'FontSize', fontsize)

else 
 
% Set title to just have end time and speed if found    
    
title(['Space $$x$$ vs. Biomass $$b$$ at $$t =$$ ' num2str(j*dt) ...
        ', $$s =$$ ' num2str(s)], 'FontSize', fontsize)
subplot(2, 1, 1)
title(['Space $$x$$ vs. Water $$w$$ at $$t =$$ ' num2str(j*dt) ...
        ', $$s =$$ ' num2str(s)], 'FontSize', fontsize) 
 
end    
    
