% Solve a particular vegetative system through continuation methods

% Zachary Singer, University of Minnesota Twin Cities, 7/28/17

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
%   Also can use data from veg_dir_sim_h1.m to compare curves.

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

% Mode Input
m = input(['Enter a mode (1 = debug, 2 = single comp, ' ...
            '3 = init. cont, 4 = tang. cont, 5 = w+ vs. s): ']);

%%% FILENAME TO SAVE DATA TO %%%
filename_std = 'st_h1m_81717_std_01.txt';
filename_bv = 'st_h1m_81717_bv_01.txt';

turn_around = 1; % 1 = on, 0 = off
turnt = 0;
        
% Step Size Constants
L = 25;             % conservative placeholder
dx = 0.05;           % change in space
ds = .5;            % change in direction
discrm_tol = .002;     % tolerance for determinant (small, > 0)
dir_num = 1000;      % number of direction searches
dir_switch = 0.25;   % switch direction after what percentage
turnt = 0;          % whether direction has switched or not
ds_turnt = 0;       % whether we have slowed down our search arclength step
x = (-L: dx: L)';   % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2);       % how many time steps

ts_data = (0 : 0.05 : 4);

% Set up Derivative and Averaging Matrix (Finite Differences)

e = ones(N, 1);
D = spdiags([-e e], 0: 1, N-1, N) / dx; % first order derivative, upwind

M = spdiags([e e], 0: 1, N-1, N) / 2;   % averaging matrix

% Initial theta and s parameter values

% Some initial values for dx / L / theta / s that work : 
% dx = 0.2, L = 25, theta = 3, s = 1
% dx = 0.1, L = 25, theta = 2.9, s = 1
% dx = 0.1, L = 25, theta = 1.5, s = 0.5
% dx = 0.05, L = 25, theta = 1.8, s = 0.5

theta = 1.8;
s = 0.5;

% Phase Condition

cnst = 5;   % steepness

b_amp = theta/(2*s) + sqrt(theta^2/(4*s^2) - 1/s); % eigenval

b_old = b_amp ./ (1+exp(cnst*x));
b_oldx =  -b_amp*cnst * exp(cnst*x) ./ (exp(cnst*x) + 1).^2; 

u0 = [b_old; b_oldx; s];

% Kinetics (uses cont_df.m)

F = @(u) cont_df(u, N, dx, theta, b_old, b_oldx, M, D);

options = optimset('Jacobian', 'on', 'Display', 'iter', 'TolFun',1e-8, ...
          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

% options = optimset('Jacobian', 'on','Display','iter','TolFun',1e-8, ...
%          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING/TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch m % mode
    
    case 1

    ord = input('Enter an order of magnitude (negative): ');
    disp(['Testing Jacobian with order 10^' num2str(ord) '...'])
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
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with s = ' num2str(s) ', theta = ' ...
            num2str(theta)])
    
       
%%%%%%%%%%%%%%%%%%% BABY CONTINUATION IN THETA: %%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    case 3

    u_new = fsolve(F, u0, options);

    b = u_new(1: N);
    v = u_new(N+1: 2*N);
    s = u_new(2*N+1);
    
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
    F = @(u) cont_df(u, N, dx, new_theta, b_old, b_oldx, M, D);
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
    
    z_new = [u_new; theta];
    
    % Set up for initial cont_dz 
    
    z_old = z_new; % initial z_old?
    dz_old = zeros(1, 2*N+2);
    dz_old(2*N+2) = 1; % initial direction in theta
    
    F = @(z) cont_dz(z, z_old, N, dx, ds, dz_old, b_old, b_oldx, M, D);
    
    % run initial cont_dz fsolve
    [z_new, fval, exitflag, output, jacobian] = fsolve(F, z_old, options);
    
    b = z_new(1: N);
    v = z_new(N+1: 2*N);
    s = z_new(2*N+1);
    
    discrm = theta^2 / (4*s^2) - 1 / s;
    
    theta = z_new(2*N+2);
    theta_data = [theta];
    s_data = [s];
    discrm_data = [discrm];
    
    z_old = z_new;
    b_old = b;
    b_oldx = v-s*b;
    
    % Update dz 
 
    Z = @(dz_var) dz_solver(dz_var, dz_old, jacobian(1:2*N+1, :), N);

    dz_new = fsolve(Z, dz_old);

    dz_new = dz_new / norm(dz_new) * ds; % normalize
    
    subplot(2, 1, 1)
    plot(x, b, 'g') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space x')
    ylabel('Biomass b')
    title(['Space x vs. Biomass b with s = ' num2str(s) ', theta = ' ...
            num2str(theta)])
    
    subplot(2, 1, 2)
    plot(s_data, theta_data, 'g')
    xlabel('Front speed s')
    ylabel('Theta')
    title(['Front speed s vs. Theta with theta = ' num2str(theta) ...
             ', s = ' num2str(s)])
             
        
    drawnow
    
    % pause
    
    plot(ts_data, 0.25*ts_data.^2, 'b')
    drawnow
    hold on
    
    cm = colormap(jet(dir_num));
    
    for i = 1: dir_num
    
    F = @(z) cont_dz(z, z_old, N, dx, ds, dz_old, b_old, b_oldx, M, D);
        
    [z_new, fval, exitflag, output, jacobian] = fsolve(F, z_old, options);
    
    b = z_new(1: N);
    v = z_new(N+1: 2*N);
    s = z_new(2*N+1);
    theta = z_new(2*N+2);
    
    discrm = theta^2 / (4*s^2) - 1 / s;
    
    
    % if getting close to 0, smaller arclength size
    % if (discrm < discrm_tol) && (turnt == 1) && (ds_turnt == 0)
    %    disp(["We are in there with discriminant = " num2str(discrm)]);
    %    ds = ds / 20;
    %    ds_turnt = 1;
    % end
    
    % stop if imaginary
    if isreal(theta) == 0 
        disp(['Theta (complex) = ' num2str(theta)])
        disp(['s (complex) = ' num2str(s)])
        break
    end   
    
    theta_data = [theta_data theta];
    s_data = [s_data s];
    discrm_data = [discrm_data discrm];
    
    b_old = b;
    b_oldx = v-s*b;
    
    if discrm < discrm_tol
        break
    end
    
    % Update dz 
 
    Z = @(dz_var) dz_solver(dz_var, dz_old, jacobian(1:2*N+1, :), N);

    dz_new = fsolve(Z, dz_old);
    
    dz_new = dz_new / norm(dz_new) * ds; % normalize, this is our dz
    
    % switch direction after a certain point (initialized above)
    
    if turn_around == 1 && i == floor(dir_switch * dir_num)
        dz_new = -1*dz_new;
        turnt = 1;
    end
           
    
    z_old = z_new + dz_new';
    dz_old = dz_new;
    
    b_height = b(floor(N/4)); % how high is interface currently
    
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
    
    %%%%%%%%%%%%%% ANALYSIS OF theta(s) %%%%%%%%%%%%%%%
      
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
    
    % Note: Generally simple is good for small number of c_values
    % you can include more colors in c_small if you want, but it's hard(er)
    % to find a good number of colors
    
    % Load text files containing data and transform continuation data    
        
    % some good files so far are:
    % 'st_h1m_81717_std_01.txt'
    cont_file = 'st_h1m_81717_std_01.txt';
    data_ct = textread(cont_file, '', 'delimiter', ',', 'emptyvalue', NaN);
    s_res = data_ct(1, :); % s data
    t_res = data_ct(2, :); % theta data

    s_res_size = size(s_res);
    s_res_size = s_res_size(2);
    
    c_results = 'c_results_82317_01.txt';       % save c range
    w_results = 'w_results_82317_01.txt';       % save w range
    mat_results = 'mat_results_82317_01.txt';   % save matrix of s values
    
    w_start = 1.1;
    d_w = 0.0001;
    w_stop = 4.5;
    
    w_range = w_start : d_w : w_stop;
    w_size = size(w_range);
    w_size = w_size(2);
    
    %%%%%%%%%%%%%%%%%%%% TRANSFORMATION AND CONSTANTS %%%%%%%%%%%%%%%%%%%%%
    
    % c parameter values (range)
    
    c = 1;     
    c_start = c;
    d_c = 1;
    % c_len = 5;  % number of c_values
    % c_stop = c + c_len*d_c - d_c;
    c_stop = 5;
    c_range = c_start : d_c : c_stop;
    
    c_new = c;
    
    % which saddle values on curve theta^2 = 4s
    s_saddle = 0 : 0.01 : 0.38;
    
    w_plus = @(c) t_res ./ sqrt(s_res + c);
    w_plus_saddle = @(c) sqrt(s_saddle ./ (s_saddle + c));
    
    % annoying setup of colors
    c_size = size(c_range);
    dir_num = c_size(2);
    cm = colormap(jet(dir_num));
    c_small = ['c' 'b' 'g' 'r' 'm'];
    c_small_num = size(c_small);
    c_small_size = c_small_num(2);
      
    
    % Where to store data
    
    res_mat = zeros(w_size, dir_num);
    eps = 1e-6; % tolerance
    
    for i = 1 : dir_num 
        
       w_res = w_plus(c_new);
       
       %%% DATA PRESERVATAION (ANNOYING) %%%
       
       s_ref = NaN(w_size, 1); % Temporary vector of NaNs
       
       % fill s_ref with w_plus values that exist
       for j = 1 : s_res_size
           w_val = w_res(j);
           s_val = s_res(j);
           idx = find(abs(w_range - w_val) < eps); % floating points
           s_ref(idx) = s_val;
       end
       
       % interpolate rest
       s_ref = fillmissing(s_ref, 'linear');
       
       res_mat(:, i) = s_ref;
        
       
       %%% OTHER PLOTTING %%%
       
       % Plotting style
       if plot_type == 1
            if i > c_small_size
                i = mod(i, 5) + 1; % stay in range of c_small 
            end
            plot(w_res, s_res, c_small(i), 'LineWidth', 2);
       
       elseif plot_type == 2
            plot(w_res, s_res, 'Color', cm(i, :)); % rainbow colors
       
       else
            warning('Must be 1 for simple or 2 for rainbow');
            
       end
       
       axis([0 6 0 6])
       xt = get(gca, 'XTick');
       set(gca, 'FontSize', tick_fontsize)
       
       hold on
       drawnow
       
       if saddle == 1
       
       w_saddle_res = w_plus_saddle(c_new);
           
       if plot_type == 1
            if i > c_small_size
                i = mod(i, 5) + 1; % stay in range of c_small 
            end
            plot(w_saddle_res, s_saddle, c_small(i), 'LineWidth', 2);
       
       elseif plot_type == 2
            plot(w_saddle_res, s_saddle, 'Color', cm(i, :)); % rainbow colors
       
       else
           warning('Must be 1 for simple or 2 for rainbow');
       
       end    
           
       ylabel('Front speed $$s$$', 'FontSize', 48)
       xlabel('Water $$w_+$$', 'FontSize', 48)
       title(['Water $$w_+$$ vs. Front speed $$s$$'], 'FontSize', 48)
       
       end
       
       c_new = c_new + d_c;
       
       hold on
       drawnow
       
    end
    
    % Scatter plots of data from direct simulation (veg_dir_sim_h1.m)
    % Commented values are additional (off curve) data points
    
    % The color coordination will work if c values match up
    
    % c = 1
    scatter([2 3 4], [1.46 2.52 3.478], c_small(1), 'LineWidth', 4);
    % 5 6, 4.38 5.44
    % c = 2
    scatter([1 2 3 4], [.584 1.82 2.857 3.81], c_small(2), 'LineWidth', 4); 
    % 5 6, 4.776 5.614
    % c = 3
    scatter([1 2 3 4], [.87 2.065 3.125 4.098], c_small(3), 'LineWidth', 4); 
    % 5 6, 5 5.926
    % c = 4
    scatter([1 2 3 4], [1.02 2.26 3.33 4.324], c_small(4), 'LineWidth', 4); 
    % 5 6, 5.246 6.25
    % c = 5
    scatter([1 2 3 4], [1.15 2.449 3.517 4.571], c_small(5), 'LineWidth', 4);  
    
    
    %%%%%%%%%%%%%%%%%%%% DATA SAVING FOR FUTURE USE %%%%%%%%%%%%%%%%%%%%%%%
    
    
    dlmwrite(c_results, c_range);
    dlmwrite(w_results, w_range);
    dlmwrite(mat_results, res_mat);
    
    otherwise
        
    warning(['Unexpected mode: Type 1 for debug, ' ...
              '2 for initial newton, 3 for baby continuation, ' ...
              '4 for tangent continuation, 5 for data transformations']);
        
end
