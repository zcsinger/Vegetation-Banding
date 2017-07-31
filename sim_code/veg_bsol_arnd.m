% Solve a particular vegetative system through continuation methods

% Zachary Singer, University of Minnesota Twin Cities, 7/20/17

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 = b_{xx} + s*b_x - m*b + b^2*(theta - s/(c+s)*b) - tau/(c+s)b_x*b^2 

% Parameters: 
% b - biomass density and function we are solving for
% m - biomass life/death constant
% db - diffusivity of biomass (scaled to 1)

% c - advective parameter
% s - wavespeed
% tau - dummy variable for continuation to 1 (originally s)
% theta - constant of integration

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
set(gcf, 'Position', [100, 800, 800, 800])
plottrue = 1;
if plottrue
    h = figure(1);
    drawnow
end

%%% Constants in Equation %%%
m = 1;    % unsure of the true range, to be safe
c = 1;      % placeholders
s = 0;
tau = 1;
%  theta = maxwell(m, c, s) % satisfy maxwell condition
theta = 1.38;

L = 50;     % conservative placeholder
dx = 0.1;
x = (-L: dx: L)'; % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2); % how many time steps

%%%%%%%%%%%%%%%%%%% Discretizing Derivatives: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = ones(N, 1);

% upwinding (forward time) for derivative D
D = spdiags([-e e], 0:1, N, N) / dx / 2;

% Neumann boundary conditions
D(1, 1) = -1/dx;
D(N, N) = 0/dx;

% central space for second derivative D2
D2 = spdiags([e -2*e e], -1:1, N, N) / dx^2;

% Neumann boundary conditions
D2(1, 1) = -1/dx^2;
D2(N, N) = -1/dx^2;

% initial point for fsolve

x0 = N/2 * dx;          % midpoint
z_init = ones(N, 1);    % artificially create heaviside


%  z_init(1 : x0 - 1) = 0;
%  z_init(x0 : N) = 1;
%  z0 = [z_init; 0]; % [b s]
%  b_old = @(x) 1 ./ (1+exp(x));           % want b_oldx (?) to be bump

%  b_oldx = @(x) exp(x) ./ (exp(x) + 1).^2; % possibly problematic

b_old = 1 ./ (1+exp(x));
b_oldx =  exp(x) ./ (exp(x) + 1).^2; 

z0 = [b_old; s];
%%%%%%%%%%%%%%%%%%% Start Iterations: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% veg_fb.m contains function F (first N comp are model, last is speed s)
F = @(z) veg_fb_arnd(z, b_old, b_oldx, D, D2, N, m, c, tau, theta, dx);
DF = @(z) veg_dfb_arnd(z, b_old, b_oldx, D, D2, N, m, c, tau, theta, dx);

% figure this out for options/diagnostics
options = optimset('Jacobian', 'off','Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',20,'Algorithm','trust-region-reflective');

% initial calculation, use this to start
z_new = fsolve(F, z0,options);
b = z_new(1:N);
s = z_new(N+1); % speed s

%%%%%%%%%%%%%%% CONTINUATION OF PARAMETERS AND PLOTTING %%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Set up plot constants for theta, c, and tau %%%%%%%%%%

% Set up theta constants
theta_start = theta;
d_theta = 0.05;
theta_stop = theta + 1.0;
theta_range = theta_start : d_theta : theta_stop;
theta_vals = NaN(size(theta_range)); % data to be filled with speeds
theta_ctr = 1; % where in theta range are we

c_start = c;
d_c = 0.05;
c_stop = c + 1.0;
c_range = c_start : d_c : c_stop;
c_vals = NaN(size(c_range)); % data to be filled with speeds
c_ctr = 1; % where in c range are we

tau_start = tau;
d_tau = 0.05;
tau_stop = tau + 1.0;
tau_range = tau_start : d_tau : tau_stop;
tau_vals = NaN(size(tau_range)); % data to be filled with speeds
tau_ctr = 1; % where in tau range are we

% PLOT INITIAL VALUES

% Plot Space x vs. Biomass b
% subplot(4, 1, 1)
plot(x, b, 'b') % space x vs biomass b
% axis([0 L 0 1])
% xlabel('Space x')
% ylabel('Biomass b')
% title(['Space x vs. Biomass b with c = ' num2str(c) ', theta = ' ...
%         num2str(theta) ', tau = ' num2str(tau)])
% % Plot theta vs. speed s
% subplot(4, 1, 2)
% plot(theta_range, theta_vals)
% % axis([theta_start theta_stop -0.1 0.05])
% xlabel('Theta')
% ylabel('Speed s')
% title('Theta vs. Speed s')
% 
% % plot c vs. s (not yet but for structure)
% subplot(4, 1, 3)
% plot(c_range, c_vals)
% % axis([c_start c_stop -0.1 0.05])
% xlabel('Advective constant c')
% ylabel('Speed s')
% title('Advective constant c vs. Speed s')
% 
% subplot(4, 1, 4)
% plot(tau_range, tau_vals)
% % axis([c_start c_stop -0.1 0.05])
% xlabel('Tau')
% ylabel('Speed s')
% title('Tau vs. Speed s')
% 
% drawnow

% pause

%%%%% Continutation Starts Below %%%%%

% for new_tau = tau_start : d_tau : tau_stop
%     F = @(z) veg_fb_arnd(z, b_old, b_oldx, D, D2, N, m, c, new_tau, theta, dx);
%     z_new = fsolve(F, z_new);
%     
%     b = z_new(1:N);
%     s = z_new(N+1);
%     tau_vals(tau_ctr) = s;
%     tau_ctr = tau_ctr + 1;
%     
%     %%%%% Plotting %%%%%
%     
%     if plottrue
%     
%     % Plot Space x vs. Biomass b
%     subplot(4, 1, 1)
%     plot(x, b, 'b') % space x vs biomass b
%      axis([0 L 0 4])
%     xlabel('Space x')
%     ylabel('Biomass b')
%     title(['Space x vs. Biomass b with c = ' num2str(c) ', theta = ' ...
%             num2str(theta) ', tau = ' num2str(new_tau)])
%         
%     % Plot theta vs. speed s
%     subplot(4, 1, 2)
%     plot(theta_range, theta_vals)
%     % axis([theta_start theta_stop -0.1 0.05])
%     xlabel('Theta')
%     ylabel('Speed s')
%     title('Theta vs. Speed s')
%     
%     subplot(4, 1, 3)
%     plot(c_range, c_vals)
%     % axis([c_start c_stop -0.1 0.05])
%     xlabel('Advective constant c')
%     ylabel('Speed s')
%     title('Advective constant c vs. Speed s')
%     
%     subplot(4, 1, 4)
%     plot(tau_range, tau_vals)
%     % axis([c_start c_stop -0.1 0.05])
%     xlabel('Tau')
%     ylabel('Speed s')
%     title('Tau vs. Speed s')
%     
%     drawnow
%      
%     % pause
%     
%     end
%     
%     %%%%% End Plotting %%%%%
%     
% end
% 
% tau = tau_stop;

%  disp('Continuing theta:')

% for new_theta = theta_start : d_theta : theta_stop
%     
%     % Change F based on new_theta value
%     F = @(z) veg_fb_arnd(z, b_old, b_oldx, D, D2, N, m, c, tau, new_theta, dx);
%     z_new = fsolve(F, z_new);
%     
%     b = z_new(1:N);
% %     b_old = b;      % possibly?
% %     b_oldx = D*b;
%     s = z_new(N+1);
%     theta_vals(theta_ctr) = s;
%     theta_ctr = theta_ctr + 1;
%     
%     %%%%% Plotting %%%%%
%     
%     if plottrue
%     
%     % Plot Space x vs. Biomass b
%     subplot(4, 1, 1)
%     plot(x, b, 'b') % space x vs biomass b
%      axis([0 L 0 4])
%     xlabel('Space x')
%     ylabel('Biomass b')
%     title(['Space x vs. Biomass b with c = ' num2str(c) ', theta = ' ...
%             num2str(new_theta) ', tau = ' num2str(tau)])
%     % Plot theta vs. speed s
%     subplot(4, 1, 2)
%     plot(theta_range, theta_vals)
%     % axis([theta_start theta_stop -0.1 0.05])
%     xlabel('Theta')
%     ylabel('Speed s')
%     title('Theta vs. Speed s')
%     
%     % plot c vs. s (not yet but for structure)
%     subplot(4, 1, 3)
%     plot(c_range, c_vals)
%     % axis([c_start c_stop -0.1 0.05])
%     xlabel('Advective constant c')
%     ylabel('Speed s')
%     title('Advective constant c vs. Speed s')
%     
%     subplot(4, 1, 4)
%     plot(tau_range, tau_vals)
%     % axis([c_start c_stop -0.1 0.05])
%     xlabel('Tau')
%     ylabel('Speed s')
%     title('Tau vs. Speed s')
%     
%     drawnow
%      
%     % pause
%     
%     end
%     
%     %%%%% End Plotting %%%%%
%                
% end
% 
% theta = theta_stop;
% disp('Continuing c:')
% 
% for new_c = c_start : d_c : c_stop
%     F = @(z) veg_fb_arnd(z, b_old, b_oldx, D, D2, N, m, new_c, tau, theta, dx);
%     z_new = fsolve(F, z_new);
%     
%     b = z_new(1:N);
%     s = z_new(N+1);
%     c_vals(c_ctr) = s;
%     c_ctr = c_ctr + 1;
%     
%     %%%%% Plotting %%%%%
%     
%     if plottrue
%     
%     % Plot Space x vs. Biomass b
%     subplot(4, 1, 1)
%     plot(x, b, 'b') % space x vs biomass b
%      axis([0 L 0 4])
%     xlabel('Space x')
%     ylabel('Biomass b')
%     title(['Space x vs. Biomass b with c = ' num2str(new_c) ', theta = ' ...
%             num2str(theta) ', tau = ' num2str(tau)])
%         
%     % Plot theta vs. speed s
%     subplot(4, 1, 2)
%     plot(theta_range, theta_vals)
%     % axis([theta_start theta_stop -0.1 0.05])
%     xlabel('Theta')
%     ylabel('Speed s')
%     title('Theta vs. Speed s')
%     
%     subplot(4, 1, 3)
%     plot(c_range, c_vals)
%     % axis([c_start c_stop -0.1 0.05])
%     xlabel('Advective constant c')
%     ylabel('Speed s')
%     title('Advective constant c vs. Speed s')
%     
%     subplot(4, 1, 4)
%     plot(tau_range, tau_vals)
%     % axis([c_start c_stop -0.1 0.05])
%     xlabel('Tau')
%     ylabel('Speed s')
%     title('Tau vs. Speed s')
%     
%     drawnow
%      
%     % pause
%     
%     end
%     
%     %%%%% End Plotting %%%%%
%     
% end
% 
% c = c_stop;

