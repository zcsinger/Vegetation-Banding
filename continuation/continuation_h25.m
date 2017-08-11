% Solve a particular vegetative system through continuation methods

% Zachary Singer, University of Minnesota Twin Cities, 8/10/17

% Take initial point from continuation and find point where theta^2 = 4s

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scaled system: 
% b_x = v - s*b
% v_x = -(theta - v)*b^2 + b

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

% Mode Input
m = input('Enter a mode (1 = debug, 2 = single comp): ');

% Step Size Constants
L = 50;             % conservative placeholder
dx = 0.05;           % change in space
x = (-L: dx: L)';   % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2);       % how many time steps

% Take theta, s values from continuation in continuation_h2.m

theta = 2.1029;
s = 1.0971;

% Set up Derivative and Averaging Matrix (Finite Differences)

e = ones(N, 1);
D = spdiags([-e e], 0: 1, N-1, N) / dx; % first order derivative, upwind
M = spdiags([e e], 0: 1, N-1, N) / 2;   % averaging matrix

cont_file = 's_theta_81017_h2_06.txt';
data_ct = textread(cont_file, '', 'delimiter', ',', 'emptyvalue', NaN);
b = data_ct(1:N, :);     % b data
v = data_ct(N+1:2*N, :); % v data

u0 = [b; v; s; theta];

% Kinetics (uses cont_df.m)

F = @(u) cont_df_h25(u, N, dx, M, D);

options = optimset('Jacobian', 'off', 'Display', 'iter', 'TolFun', 1e-8, ...
          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

% reference if you want to turn Display off (or other)
% options = optimset('Jacobian', 'on','Display','iter','TolFun',1e-8, ...
%          'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');

%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING/TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch m % mode
    
    case 1

    ord = input('Enter an order of magnitude (negative): ');
    disp(['Testing Jacobian with order 10^' num2str(ord) '...'])
    du = 10^ord * rand(2*N + 2, 1); % random points, order 'ord'
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
    theta = u_new(2*N+2);

    plot(x, b, 'b') % space x vs biomass b
    % axis([0 L 0 1])
    xlabel('Space $$x$$', 'FontWeight', 'bold', 'FontSize', 24)
    ylabel('Biomass $$b$$', 'FontWeight', 'bold', 'FontSize', 24)
    title(['Space $$x$$ vs. Biomass $$b$$ with $$\theta$$ = ' num2str(theta) ...
        ', $$s$$ = ' num2str(s)], 'FontWeight', 'bold', 'FontSize', 20)
               
    otherwise
        
    warning('Unexpected mode: Type 1 for debug, 2 for fsolve');
        
end