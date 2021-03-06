% Solve a particular vegetative system through continuation methods
% Note: Need intersections.m for mode 6

% Zachary Singer, University of Minnesota Twin Cities, 8/9/17

% Differences between continuation_h1.m and continuation_h2.m : 
%
% - Looking for lower edge heteroclinic solution. Boundary conditions change.
% - Format is (previous boundary condition), current boundary condition :
%   - (b(N) - lamb_neg*v(N)), b(1) - lamb_pos*v(1)
%   - (dot([b(1) v(1)] - [b_plus s*b_plus], stab_eval)), 
%       dot([b(N) v(N)] - [b_minus s*b_minus], unstab_eval)
%
% - Different initial interfaces 

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scaled system: 
% b_x = v - s*b
% v_x = -(theta - v)*b^2 + b

% theta - varied parameter (want to find theta in terms of s)

%%%%%%%%%%%%%%%%%%%% Specifications and Guide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Created to use tangent arclength continuation to find a more accurate
% relationship between front speed s and parameter theta in model above.

% Relies on cont_df.m for function and jacobian calculations, and
% cont_dz.m for continuation method calculations

% How to use:

% There are modes you can choose for particular tasks:
% plotting  - see plots of space x vs. biomass b and space x vs. water w
%

% Mode 1 (Debug):
%   Debugs by testing to make sure that the Jacobian expression is correct.
%   Tests different magnitudes to show quadratic error.

% Mode 2 (Initial Newton):
%   A single newton iteration to find the first solution to the nonlinear
%   equation given in cont_df.m.
%   F takes in u = [b; v; s] (2*N+1 variables). Has (N-1) + (N-1) eqns.
%   related to b and v, as well as 2 boundary cond. and 1 phase cond.
%   Plots Space x vs. Biomass b just to see interface, not too instructive. 

% Mode 3 ("Baby Continuation" in Theta):
%   Sets up initial newton iteration as in Mode 2, but linearly varies
%   Theta in one direction and plots the interface as well as 
%   Theta vs. Front Speed s 

% Mode 4 (Tangent Arclength Continuation):
%   First will find single solution using fsolve, and then will use 
%   continuation in theta with initial direction dz_old (line ~225).
%   Then will use tangent arclength continuation to find new directions
%   and new solutions. Plots space x vs. biomass b (not too useful) and
%   theta vs. front speed s. Writes data to file.

% Mode 5 (Data Transformation and Plotting w_plus vs s):
%   Given data found in Mode 4, transforms the data by using scalings
%   from our model back to modified Klausmeier, and plots the curves for
%   the height of the water front w_plus vs. front speed s.
%   Also can use data from veg_model_80317.m to compare curves.

% Mode 6 (Bifurcation Diagram for theta - s plane): 
%   You can create diagrams with any combination of the saddle curve, 
%   Hopf curve, and both heteroclinic curves before/on the saddle curve.

% Mode 7 (Heteroclinic Intersection) : 
%   Transforms the lower edge heteroclinic to plot w_- vs. s_lower instead of 
%   b_+ vs. s. On the same plot is the relationship between showing
%   w_+ vs. s^upper and w_- vs. s_lower. Dashed arrow added manually.

% Suggested Ranges of Variables:
% dx : 0.1 to 0.5 (smaller better)
% theta / s : a list of converging initial conditions given (picky)
% b_old / b_oldx : phase condition, looks like step function 
% dz_old (in mode 4) : initial direction solely in direction of theta
%                      (i.e. (0, 0, ..., 0, 1))

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

tick_fontsize = 36;
fontsize = 48;
linewidth = 4;

% Mode Input
m = input(['Enter mode (Press 8 for help): ']);

% Options for turning around or slowing down (turn "on" or "off")
        
slow_down = "off";  % turn on for slowdown near 0 discriminant
discrm_tol = .05;     % tolerance for determinant (small, > 0)
ds_turnt = 0;       % whether we have slowed down our search arclength step
ds_slowdown = 0.2;  % percentage of ds after slowdown

turn_around = "on"; % turn on for turnaround after dir_switch % of dir_num
dir_switch = 0.05;   % switch direction after what percentage
turnt = 0;          % whether direction has switched or not


% !!! FILENAME TO SAVE DATA TO !!!
filename_std = 's_theta_date_std_h2m_01.txt';
filename_bv = 's_theta_date_bv_h2m_01.txt';

        
% Step Size Constants
L = 25;             % conservative placeholder
dx = 0.05;           % change in space
ds = .1;            % change in direction
dir_num = 500;      % number of direction searches
x = (-L: dx: L)';   % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2);       % how many time steps

ts_data = (0 : 0.05 : 3);

% Initial theta and s parameter values

% Some initial values for dx / L / theta / s that work : 
% dx = 0.2, L = 25, theta = 1.7, s = 0.5
% dx = 0.05, L = 25, theta = 1.8, s = 0.5
% dx = 0.05, L = 25, theta = 1.5, s = 0.4

theta = 1.5;
s = 0.4;


% Set up Derivative and Averaging Matrix (Finite Differences)

e = ones(N, 1);
D = spdiags([-e e], 0: 1, N-1, N) / dx; % first order derivative, upwind

M = spdiags([e e], 0: 1, N-1, N) / 2;   % averaging matrix

% Phase Condition

cnst = 5;   % steepness

b_amp = theta/(2*s) + sqrt(theta^2/(4*s^2) - 1/s); % eigenval (CHANGE??)

% b_old = -b_amp ./ (1+exp(cnst*x));
b_old = (b_amp / 2) * tanh(cnst*x) + (b_amp / 2);
% b_oldx =  b_amp*cnst * exp(cnst*x) ./ (exp(cnst*x) + 1).^2;
b_oldx = (b_amp / 2) * cnst * (sech(cnst*x)).^2;

u0 = [b_old; b_oldx; s];

% initial conditions for Jacobian (z) test

z_old = [b_old; b_oldx; s; theta];
z0 = z_old;
dz_old = zeros(1, 2*N+2);
dz_old(2*N+2) = 1; % initial direction in theta

% Kinetics (uses cont_df.m)

F = @(u) cont_df_h2(u, N, dx, theta, b_old, b_oldx, M, D);

FZ = @(z) cont_dz_h2(z, z_old, N, dx, ds, dz_old, b_old, b_oldx, M, D);

options = optimset('Jacobian', 'on', 'Display', 'iter', 'TolFun', 1e-8, ...
          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

% reference if you want to turn Display off (or other)
% options = optimset('Jacobian', 'on','Display','iter','TolFun',1e-8, ...
%          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING/TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch m % mode
    
    case 0 % test jacobian with b, v, s, theta and dz
        
    ord = input('Enter an order of magnitude (negative): ');
    disp(['Testing Jacobian (z) with order 10^' num2str(ord) '...'])
    dz = 10^ord * rand(2*N + 2, 1); % random points, order 'ord'
    [F0, J0] = FZ(z0);
    [F1, J1] = FZ(z0 + dz);
    err_z = (F1 - F0 - J0*dz);
    disp(['Error: ' num2str(norm(err_z))]);
    plot(err_z)
    
    case 1 % test initial jacobian with b, v, s

    ord = input('Enter an order of magnitude (negative): ');
    disp(['Testing Jacobian (u) with order 10^' num2str(ord) '...'])
    du = 10^ord * rand(2*N + 1, 1); % random points, order 'ord'
    [F0, J0] = F(u0);
    [F1, J1] = F(u0 + du);
    err_u = (F1 - F0 - J0*du);
    disp(['Error: ' num2str(norm(err_u))]);
    plot(err_u)


%%%%%%%%%%%%%%%%%%% SINGLE USE OF FSOLVE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2

    u_new = fsolve(F, u0, options);

    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);

    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space $$x$$', 'FontSize', 24)
    ylabel('Biomass $$b$$', 'FontSize', 24)
    title(['Space $$x$$ vs. Biomass $$b$$ with $$s =$$ ' num2str(s) ', $$\theta =$$ ' ...
            num2str(theta)], 'FontSize', 20)
    
       
%%%%%%%%%%%%%%%%%%% BABY CONTINUATION IN THETA: %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    case 3

    u_new = fsolve(F, u0, options);

    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);
    
    % Initial Plot
    
    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space $$x$$', 'FontSize', 24)
    ylabel('Biomass $$b$$', 'FontSize', 24)
    title(['Space $$x$$ vs. Biomass $$b$$ with $$s =$$ ' num2str(s) ', $$\theta_0 =$$ ' ...
            num2str(theta)], 'FontSize', 20)
    
    pause    
        
    % update phase condition
    
    b_old = b;
    b_oldx = v - s*b;
    
    % set up theta arrays for plotting later
    
    theta_start = theta;
    d_theta = 0.05;
    theta_stop = theta - 1.5;
    % set theta_go/stop so direction of continuation d/n matter
    if theta_start > theta_stop
        d_theta = -d_theta;
    end
    theta_go = min(theta_start, theta_stop);
    theta_end = max(theta_start, theta_stop);
    % set range of theta for plotting
    theta_range = theta_start : d_theta : theta_stop;
    theta_vals = NaN(size(theta_range)); % data to be filled with speeds
    theta_ctr = 1; % where in tau range are we

    % initial theta plotting
    
    subplot(2, 1, 2)
    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with s = ' num2str(s) ', theta = ' ...
            num2str(theta)])
    subplot(2, 1, 2)
    plot(theta_range, theta_vals)
    axis([theta_go theta_end 0 2])
    xlabel('Theta')
    ylabel('Speed s')
    title('Theta vs. Speed s')
 
    drawnow
    
    pause
        
    % vary theta
    
    for new_theta = theta_start : d_theta : theta_stop
    
    % Change F based on new_theta value
    F = @(u) cont_df_h2(u, N, dx, new_theta, b_old, b_oldx, M, D);
    u_new = fsolve(F, u_new, options);
    
    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);
    
    b_old = b;
    b_oldx = v - s*b;
    
    theta_vals(theta_ctr) = s;
    theta_ctr = theta_ctr + 1;
    
    %%%%% Plotting %%%%%
    
    % Plot Space x vs. Biomass b
    subplot(2, 1, 1)
    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 4])
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with theta = ' num2str(new_theta) ...
            ', s = ' num2str(s)])
    % Plot theta vs. speed s
    subplot(2, 1, 2)
    plot(theta_range, theta_vals)
    axis([theta_go theta_end 0 2])
    xlabel('Theta')
    ylabel('Speed s')
    title('Theta vs. Speed s')
 
    drawnow
   
    end
  
%%%%%%%%%%%%%%%%%%%%%% TANGENT ARCLENGTH CONTINUATION %%%%%%%%%%%%%%%%%%%%%

    case 4
        
    % number of direction iterations
    theta_data = [];
    s_data = []; % want to plot theta vs. s
    
        
    u_new = fsolve(F, u0, options);

    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);
    
    % Plot after u fsolve
    
    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space $$x$$', 'FontSize', 24)
    ylabel('Biomass $$b$$', 'FontSize', 24)
    title(['(After $$u$$ fsolve) Space $$x$$ vs. Biomass $$b$$ with $$s =$$ ' ...
            num2str(s) ', $$\theta_0 =$$ ' num2str(theta)], 'FontSize', 20)
    
    pause 
    
    z_new = [u_new; theta]; 
    
    % Set up for initial cont_dz 
    
    z_old = z_new; % initial z_old?
    dz_old = zeros(1, 2*N+2);
    dz_old(2*N+2) = 1; % initial direction in theta
    
    F = @(z) cont_dz_h2(z, z_old, N, dx, ds, dz_old, b_old, b_oldx, M, D);
    
    % run initial cont_dz fsolve
    [z_new, fval, exitflag, output, jacobian] = fsolve(F, z_old, options);
    
    b = z_new(1: N);
    v = z_new(N+1: 2*N);
    s = z_new(2*N+1);
    theta = z_new(2*N+2);

    theta_data = [theta];
    s_data = [s];
    
    discrm = theta^2 / (4*s^2) - 1 / s;
    discrm_data = [discrm];
    
    z_old = z_new;
    b_old = b;
    b_oldx = v-s*b;
    
    % Update dz 
 
    Z = @(dz_var) dz_solver(dz_var, dz_old, jacobian(1:2*N+1, :), N);

    dz_new = fsolve(Z, dz_old);

    dz_new = dz_new / norm(dz_new) * ds; % normalize
    
    % Initial Plot
    
    subplot(2, 1, 1)
    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space $$x$$', 'FontSize', 24)
    ylabel('Biomass $$b$$', 'FontSize', 24)
    title(['(After $$z$$ fsolve) Space $$x$$ vs. Biomass $$b$$ with $$s =$$ ' ...
            num2str(s) ', $$\theta_0 =$$ ' num2str(theta)], 'FontSize', 20)

    subplot(2, 1, 2)
    plot(s_data, theta_data, 'g')
    xlabel('Front speed $$s$$', 'FontSize', 24)
    ylabel('$$\theta$$', 'FontSize', 24)
    title(['(After $$z$$ fsolve) Front speed $$s$$ vs. $$\theta$$ with $$\theta =$$ ' ...
             num2str(theta) ', $$s =$$ ' num2str(s)], 'FontSize', 20)
            
    % curve theta^2 = 4s
    plot(ts_data, 0.25*ts_data.^2, 'b')
    drawnow
    hold on
    
    pause
    
       
    cm = colormap(jet(dir_num)); % rainbow coloring
    
    for i = 1: dir_num
    
    F = @(z) cont_dz_h2(z, z_old, N, dx, ds, dz_old, b_old, b_oldx, M, D);
        
    [z_new, fval, exitflag, output, jacobian] = fsolve(F, z_old, options);
    
    b = z_new(1: N);
    v = z_new(N+1: 2*N);
    s = z_new(2*N+1);
    theta = z_new(2*N+2);
    
    discrm = theta^2 / (4*s^2) - 1 / s; % discriminant of equilibrium
    
    disp(['discrim = ' num2str(discrm) ', theta = ' num2str(theta) ...
            ', s = ' num2str(s)]);
    
    if slow_down == "on"
        % if getting close to 0, smaller arclength size
        if (discrm < discrm_tol) && (ds_turnt == 0)
            disp(["We are in there with discriminant = " num2str(discrm)]);
            ds = ds_slowdown * ds;
            ds_turnt = 1;
        end
    end
    
    % stop if imaginary
    if isreal(theta) == 0 
        disp(['Theta (complex) = ' num2str(theta)])
        disp(['s (complex) = ' num2str(s)])
        break
    end   
    
    % Update data and variables
    
    theta_data = [theta_data theta];
    s_data = [s_data s];
    discrm_data = [discrm_data discrm];
    
    b_old = b;
    b_oldx = v-s*b;
    
    % Update dz 
 
    Z = @(dz_var) dz_solver(dz_var, dz_old, jacobian(1:2*N+1, :), N);

    dz_new = fsolve(Z, dz_old);
    
    dz_new = dz_new / norm(dz_new) * ds; % normalize, this is our dz
    
    % switch direction after a certain point (initialized above)
    
    % i == floor(dir_switch * dir_num) 
    if turn_around == 'on' && i == floor(dir_switch * dir_num)
        dz_new = -1*dz_new;
        turnt = 1;
    end
    
    z_old = z_new + dz_new';
    dz_old = dz_new;
    
    b_height = b(floor(3*N/4)); % how high is interface currently
    
    subplot(2, 1, 1)
    plot(x, b, 'Color', cm(i, :)) % space x vs biomass b
    xlabel('Space $$x$$', 'FontWeight', 'bold', 'FontSize', 24)
    ylabel('Biomass $$b$$', 'FontWeight', 'bold', 'FontSize', 24)
    title(['Space $$x$$ vs. Biomass $$b$$ with $$b_+$$ = ' ...
             num2str(b_height)], 'FontWeight', 'bold', 'FontSize', 20)
    hold on
    
    subplot(2, 1, 2)
    plot(theta_data, s_data, 'r')
    ylabel('Front speed $$s$$', 'FontWeight', 'bold', 'FontSize', 24)
    xlabel('$$\theta$$', 'FontWeight', 'bold', 'FontSize', 24)
    title(['$$\theta$$ vs. Front speed $$s$$ with $$\theta$$ = ' num2str(theta) ...
             ', $$s$$ = ' num2str(s) ', discrim = ' num2str(discrm)], ...
             'FontWeight', 'bold', 'FontSize', 20)
    % axis([1.4 2.6 0 1.5])
    hold on
    
    drawnow
    
    end
    
    %%%%%%%%%%%%%% Set up figures and write data to file %%%%%%%%%%%%%%%
    
    % change plot titles
    title('$$\theta$$ vs. Front speed $$s$$', ...
            'FontWeight', 'bold', 'FontSize', 24)
    subplot(2, 1, 1)
    title('Space $$x$$ vs. Biomass $$b$$', ...
          'FontWeight', 'bold', 'FontSize', 24)
      
    results1 = [s_data; theta_data; discrm_data];
    results2 = [b; v];
    
    dlmwrite(filename_std, results1);
    dlmwrite(filename_bv, results2);
    
%%%%%%%%%%%%%%%%%%% SCALING AND DATA TRANSFORMATIONS %%%%%%%%%%%%%%%%%%%%%%

    case 5
           
    %%%%% MODES (PLOT TYPE, INCLUDING SADDLE CURVE) %%%%%    
    
    plot_type = input('Plot Type (input 1 for simple, 2 for rainbow): ');    
    saddle = input('Saddle (input 1 for inclusion of saddle trans.): ');
    
    fontsize = input('Fontsize (i.e. 48): ');
    tick_fontsize = input('Tick Fontsize (i.e. 36): ');
    linewidth = input('LineWidth (i.e. 4): ');
    
    % Note: Generally simple is good for small number of c_values
    % you can include more colors in c_small if you want, but it's hard(er)
    % to find a good number of colors  
        
    % Load text files containing data and transform continuation data    
        
    % some good files so far are:
    % 's_theta_81417_std_h2m_01.txt'
    cont_file = 's_theta_81417_std_h2m_01.txt';
    data_ct = textread(cont_file, '', 'delimiter', ',', 'emptyvalue', NaN);
    s_res = data_ct(1, :); % s data
    t_res = data_ct(2, :); % theta data
    
%%%%%%%%%%%%%%%%%%%%%%%% TRANSFORMATION AND CONSTANTS %%%%%%%%%%%%%%%%%%%%%
    
    %%% Parameter of c values %%%
    c = 1;  
    c_start = c;
    d_c = 1;
    c_stop = 5;
    c_range = c_start : d_c : c_stop;
    
    c_new = c;
    
    %%% Where we draw saddle curve as part of heteroclinic curve %%%
    
    s_saddle = 1.25 : .025 : 3;
   
    %%% Transformations %%%
    
    w_minus = @(c) t_res ./ sqrt(s_res + c);
    
    b_plus = @(c) (sqrt((s_res + c)) ./ (2*s_res)) .* ...
              (t_res + sqrt(t_res.^2 - 4*s_res));
    
    b_plus_saddle = @(c) sqrt((s_saddle + c) ./ s_saddle);
        
          
    % annoying setup for c_values
    c_size = size(c_range);
    dirs = c_size(2);
    cm = colormap(jet(dirs));
    c_small = ['c' 'b' 'g' 'r' 'm'];
    c_small_num = size(c_small);
    c_small_size = c_small_num(2);
    
    % go through each c value desired
    for i = 1 : dirs

       b_res = b_plus(c_new);
       
       % Plotting style (simple is c_small values, rainbow is all colors)
       if plot_type == 1
            if i > c_small_size
                i = mod(i, 5) + 1; % stay in range of c_small 
            end
            plot(b_res, s_res, c_small(i), 'LineWidth', linewidth);
       
       elseif plot_type == 2
            plot(b_res, s_res, 'Color', cm(i, :), 'LineWidth', linewidth); 
            % rainbow colors
       
       else
            warning('Must be 1 for simple or 2 for rainbow');
            
       end
       
       axis([0 10 0 2])
       xt = get(gca, 'XTick');
       set(gca, 'FontSize', tick_fontsize)
       
       hold on
       drawnow
       
       if saddle == 1
       
       b_saddle_res = b_plus_saddle(c_new);
           
       if plot_type == 1
            if i > c_small_size
                i = mod(i, 5) + 1; % stay in range of c_small 
            end
            plot(b_saddle_res, s_saddle, c_small(i), 'LineWidth', linewidth);
       
       elseif plot_type == 2
            plot(b_saddle_res, s_saddle, 'Color', cm(i, :), ...
                    'LineWidth', linewidth); % rainbow colors
       
       else
           warning('Must be 1 for simple or 2 for rainbow');
       
       end    
           
       ylabel('Front speed $$s$$', 'FontSize', fontsize)
       xlabel('Biomass $$b_+$$', 'FontSize', fontsize)
       title(['Biomass $$b_+$$ vs. Front speed $$s$$'], 'FontSize', fontsize)
       
       end
       
       c_new = c_new + d_c;
       
       hold on
       drawnow
       
       
    end
    
    % Scatter plots of data from direct simulation (veg_dir_sim_h1.m)
    
    % c = 1
    scatter([2 3 4], [0.7713 0.5114 0.3264], c_small(1), ...
            'LineWidth', linewidth) % 5, 0.016
    % c = 2
    scatter([2 3 4 5], [0.95416 0.65 0.485 0.39104], c_small(2), ...
            'LineWidth', linewidth)
    % c = 3
    scatter([2 3 4 5], [1.0837 0.7587 0.5767 0.4637], c_small(3), ...
            'LineWidth', linewidth) % 2, 1.0837
    % c = 4
    scatter([2 3 4 5], [1.1830 0.8453 0.6481 0.5238], c_small(4), ...
            'LineWidth', linewidth) % 2, 1.1830
    % c = 5
    scatter([2 3 4 5], [1.362 0.9186 0.7095 0.5759], c_small(5), ...
            'LineWidth', linewidth) % 2, 1.362
    
%%%%%%%%%%%%%%%%%%%%%%%% REPLOTTING IN S-THETA PLANE %%%%%%%%%%%%%%%%%%%%%%

    
    case 6
           
    % Load text files for s and theta plot (to generate figure)    
       
    fnt = input('Font Size: ');
    tick_fontsize = input('Tick Font Size: ');
    
    % some good files so far for lower edge het2 are:
    % 's_theta_81417_std_h2m_01.txt'
    cont_file_h2 = 's_theta_81417_std_h2m_01.txt';
    data_ct_h2 = textread(cont_file_h2, '', 'delimiter', ',', 'emptyvalue', NaN);
    s_res_h2 = data_ct_h2(1, :); % s data
    t_res_h2 = data_ct_h2(2, :); % theta data
    
    % some good files so far for upper edge het1 are:
    % 'st_h1m_81617_std_01.txt'
    cont_file_h1 = 'st_h1m_81617_std_01.txt';
    data_ct_h1 = textread(cont_file_h1, '', 'delimiter', ',', 'emptyvalue', NaN);
    s_res_h1 = data_ct_h1(1, :); % s data
    t_res_h1 = data_ct_h1(2, :); % theta data
    
    % print where intersections are
    [x_intr, y_intr] = intersections(t_res_h2, s_res_h2, t_res_h1, s_res_h1, 1);
    
    ts_1 = 0 : 0.01 : 1.23;
    ts_2 = 2.2 : 0.01 : 3;
    
    hopf = 0 : 0.01 : 1;
    theta_hopf = (hopf.^2 + 1) ./ sqrt(hopf);
    
    
    % Uncomment this for the saddle curve
    plot(ts_data, 0.25*ts_data.^2, 'b', 'LineWidth', 4); % theta^2 = 4s
    
    axis([0 3 0 3])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', tick_fontsize)
    hold on
    
    % Uncomment this for the saddle curve intersection with upper edge
    % plot(ts_1, 0.25*ts_1.^2, 'g', 'LineWidth', 2);
    
    % Uncomment this for the saddle curve intersection with lower edge
    % plot(ts_2, 0.25*ts_2.^2, 'r', 'LineWidth', 2);
    
    % Uncomment this for the hopf curve, -- for dashed
    plot(theta_hopf, hopf, '--m', 'LineWidth', 4);
    
    % Uncomment this for the lower edge heteroclinic before saddle curve
    % plot(t_res_h2, s_res_h2, 'r', 'LineWidth', 2); % s vs. theta for h2
    
    % Uncomment this for the upper edge heteroclinic before saddle curve
    % plot(t_res_h1, s_res_h1, 'g', 'LineWidth', 2); % s vs. theta for h1
    
    % Uncomment this for approximate intersection point
    % scatter([1.8464], [0.7688], 'k', 'filled')
    
    % Plot Labels / Legend
    xlabel('$$\theta$$', 'FontSize', fnt)
    ylabel('Front speed $$s$$', 'FontSize', fnt)
    title('$$\theta$$ vs. Front speed $$s$$', 'FontSize', fnt)
    
    % You will need to fix the legend depending on what you decide to plot
    legend({'$$\theta^2 = 4s$$ (Saddle)', '$$\theta = \frac{s^2 + 1}{\sqrt{s}}$$ (Hopf)'}, ... 
        'Interpreter', 'latex', 'FontSize', fnt, 'Location', 'northwest')
    legend('boxoff')
    
    % Quick reference for hopf curve equation
    % '$$\theta = \frac{s^2 + 1}{\sqrt{s}}$$'}, ...
    
    
%%%%%%%%%%%%%%%%%%%%%%%% HETEROCLINIC INTERACTION %%%%%%%%%%%%%%%%%%%%%%

   
    case 7
   
    % Load text files for s and theta plot (to generate figure)    
       
    fontsize = input('Font Size: ');
    tick_fontsize = input('Tick Font Size: ');
    
    % some good files so far for lower edge het2 are:
    % 's_theta_81417_std_h2m_01.txt'
    cont_file_h2 = 's_theta_81417_std_h2m_01.txt';
    data_ct_h2 = textread(cont_file_h2, '', 'delimiter', ',', 'emptyvalue', NaN);
    s_res_h2 = data_ct_h2(1, :); % s data
    t_res_h2 = data_ct_h2(2, :); % theta data
    
    % some good files so far for upper edge het1 are:
    % 'st_h1m_81617_std_01.txt'
    cont_file_h1 = 'st_h1m_81617_std_01.txt';
    data_ct_h1 = textread(cont_file_h1, '', 'delimiter', ',', 'emptyvalue', NaN);
    s_res_h1 = data_ct_h1(1, :); % s data
    t_res_h1 = data_ct_h1(2, :); % theta data
    
    c = 1;  % first value will be c + 1
    c_start = c;
    d_c = 0.5;
    c_stop = 5;
    
    c_new = c;
   
    % initial test
    w_plus_h1 = @(c) t_res_h1 ./ sqrt(s_res_h1 + c);
    w_plus_h2 = @(c) t_res_h2 ./ sqrt(s_res_h2 + c);
    
    c_range = size(c_start : d_c : c_stop);
    dir_num = c_range(2);
    cm = colormap(jet(dir_num));
    c_small = ['c' 'b' 'g' 'r' 'm'];
         
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', tick_fontsize)
    
    hold on 
    
    for i = 1 : dir_num
          
       c_new = c_new + d_c;
       
       w_res_h1 = w_plus_h1(c_new);
       w_res_h2 = w_plus_h2(c_new);
       
       % plot(w_res_h1, s_res_h1, c_small(i), 'LineWidth', 2)
        plot(w_res_h1, s_res_h1, 'Color', cm(i, :), 'LineWidth', 2) % rainbow colors
       
       % plot(w_res_h2, s_res_h2, c_small(i), 'LineWidth', 2)
        plot(w_res_h2, s_res_h2, 'Color', cm(i, :), 'LineWidth', 2); 
       
       ylabel('Front speed $$s^u / s^l$$', 'FontSize', fontsize)
       xlabel('Water $$w_+ / w_-$$', 'FontSize', fontsize)
       title(['Water $$w_+ / w_-$$ vs. Front speed $$s^u / s^l$$'], 'FontSize', fontsize)
       
       drawnow
       
    end
    
    case 8
        
    fprintf(['0 or 1 - debug for jacobian, 2 - single computation, \n' ...
            '3 - baby continuation, 4 - tangent continuation, \n' ...
            '5 - transformations to b_+ vs. s, 6 - bifurcation diagrams, \n' ...
            '7 - heteroclinic interaction. \n']); 
    
    otherwise
        
    warning(['Unexpected mode: Type 0/1 for debug, ' ...
              '2 for init, 3 for baby cont, ' ...
              '4 for tang cont, 5 for data transf, 6 for figure, ' ...
              '7 for heteroclinic interaction.']);
        
end
