% Simulation of modified Klausmeier model of two-species diffusion-reaction
% Looking for homoclinic orbits.

% Zachary Singer, University of Minnesota - Twin Cities, 8/16/17

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


%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Original system with traveling wave ansatz: 
% b_t = d_b*b_xx + w*b^2 - mb + sb_x
% w_t =          - w*b^2 + mb + (s+c)w_x

% Convention : d_b = 1, m = 1 (for simplicity)

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
speed = "off";
speed_alt = "off";
write_data = "off";
filename = "veg_model_hom_81617_01.txt";

fontsize = 36;
tick_fontsize = 24;
linewidth = 4;

pause_time = 3.5;

%%%%%% IMPORTANT CONSTANTS %%%%%%

% get theta, s from homoclinic curve (arnd's)

theta = 2;
s = 0.25;
c = 2;
L_wid = 2;          % how spaced out mid is (2*L_wid)
tmax = 10;          % how long time runs for

% Other Constants 

w_plus = theta/sqrt(s+c);
w_orig = w_plus - 1 / w_plus;

w_plus_s = w_plus;        % starting point for b_plus
d_w = 1;
w_plus_e = w_plus; % ending point for b_plus
w_range = w_plus_s : d_w : w_plus_e;

w_data = w_range; % this could be fishy
s_data = [];

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

t_eps = 0.1;    % threshold for front when calculating speed s
t_start = 0;    % set to 0 for now
t_stop = 0;     % placeholder
t_size = floor(L/10) / dx;    % !!! Adjust this to change speed accuracy !!!
t_thres = floor(L/20) / dx;  % !!! Amount of time to wait until t_start !!!
s = 0;          % placeholder

% Derivative Matrices

e = ones(N, 1);

D1 = spdiags([-e e], 0 : 1, N, N) / dx; % first derivative, upwind
% boundary conditions? --> see veg_model_f.m F(2*N) component
D2 = spdiags([e -2*e e], -1 : 1, N, N) / dx^2; % second deriv., central
% boundary conditions?

%%%%%%%%%%%%%%%%%%% Loop over w_plus values: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for w_plus = w_range

b_plus = 1 / w_plus;
    
% Kinetics (found in veg_model_f.m and veg_model_j.m)

F = @(t, V) veg_model_fp(t, V, D1, D2, N, dx, c, w_plus);
J = @(t, V) veg_model_j(t, V, D1, D2, N, c);

options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'Jacobian', J);

% Interfaces for biomass b and water w

const = 5; % steepness of interface in middle

% Opposite interface as in veg_model_80317.m
b = ((w_plus / 2) * tanh(const*(x - floor(L/2 - L_wid)))) + w_plus / 2 - ...
     (((w_plus / 2) * tanh(const*(x - floor(L/2 + L_wid)))) + w_plus / 2);
w = (-(w_orig / 2) * (tanh(const*(x - floor(L/2 - L_wid))))) + w_orig / 2 + 1 / w_plus + ...
     (((w_orig / 2) * tanh(const*(x - floor(L/2 + L_wid)))) + w_orig / 2); 

U = [b; w];

if plotting == "on"

% Initial Plot to view Interfaces

subplot(2, 1, 1)
plot(x, w, 'b', 'LineWidth', linewidth)    
axis([0 L 0 w_plus + 4])
xt = get(gca, 'XTick');
set(gca, 'FontSize', tick_fontsize) 
xlabel('Space $$x$$', 'FontSize', fontsize)
ylabel('Water $$w$$', 'FontSize', fontsize)
title(['Initial Space $$x$$ vs. Water $$w$$ at time $$t = 0$$'], ...
        'FontSize', fontsize)

subplot(2, 1, 2)
plot(x, b, 'g', 'LineWidth', linewidth)
axis([0 L 0 w_plus + 4]) % Might have to increase y-axis
xlabel('Space $$x$$', 'FontSize', fontsize)
ylabel('Biomass $$b$$', 'FontSize', fontsize)
title(['Initial Space $$x$$ vs. Biomass $$b$$ at time $$t = 0$$'], ...
        'FontSize', fontsize)

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
        
    % speed calculation
    w_sp = find(w > w_plus - t_eps);
    t_sp = find(w_sp > floor(N/2));
    w_sp = w_sp(t_sp); % find indices greater than floor(N/2)
    wf_size = size(w_sp);
    wf_len = wf_size(1);
    
    if t_start == 0 && wf_len > 0
        if w_sp(1) > (floor(N/2) + t_thres)
            w_sp(1);
            t_start = j*dt;
        end
    elseif t_stop == 0 && wf_len > 0 
        if w_sp(1) > (floor(N/2) + t_thres + t_size) 
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
    axis([0 L 0 w_plus + 4])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', tick_fontsize) 
    xlabel('Space $$x$$', 'FontSize', fontsize)
    ylabel('Water $$w$$', 'FontSize', fontsize)
    title(['Space $$x$$ vs. Water $$w$$ at time $$t =$$ ' num2str(j*dt)], ...
            'FontSize', fontsize)
    
    subplot(2, 1, 2)
    plot(x, b, 'g', 'LineWidth', linewidth)
    axis([0 L 0 w_plus + 4])
    xlabel('Space $$x$$', 'FontSize', fontsize)
    ylabel('Biomass $$b$$', 'FontSize', fontsize)
    title(['Space $$x$$ vs. Biomass $$b$$ at time $$t =$$ ' ...
             num2str(j*dt)], 'FontSize', fontsize)
        
    drawnow
    % pause
    
    if j*dt == pause_time
        pause
    end
    
    end
    
end

if speed_alt == 'on'
   w = U(N+1 : 2*N);
   w_minus = w(floor(N/5)); % look at flat section
   w_plus = 1 / b_plus;
   
   s = c * (w_minus - w_plus) / (b_plus + w_plus - w_minus);
   s_data = [s_data s];
end

title(['Space $$x$$ vs. Biomass $$b$$ at $$t =$$ ' num2str(j*dt)], ...
        'FontSize', fontsize)
subplot(2, 1, 1)
title(['Space $$x$$ vs. Water $$w$$ at $$t = $$ ' num2str(j*dt)], ...
        'FontSize', fontsize)


end

%%%%%% WRITE DATA (turn write_data = "on" above) %%%%%

if write_data == "on"

results = [w_data; s_data];
dlmwrite(filename, results);

end
