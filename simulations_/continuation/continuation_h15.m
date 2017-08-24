% Solve a particular vegetative system through continuation methods

% Zachary Singer, University of Minnesota Twin Cities, 8/16/17

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

N = 158; % CHANGE THIS TO MATCH DATA
dx = 0.05;

% Good previous values: theta = 1.3477, s = 0.45367

theta = 1.3477;
s = 0.45367;

% Set up Derivative and Averaging Matrix (Finite Differences)

e = ones(N, 1);
D = spdiags([-e e], 0: 1, N-1, N) / dx; % first order derivative, upwind
M = spdiags([e e], 0: 1, N-1, N) / 2;   % averaging matrix

cont_file = 'st_h1m_81617_bv_01.txt';
data_ct = textread(cont_file, '', 'delimiter', ',', 'emptyvalue', NaN);
b = data_ct(1:N, :);     % b data
v = data_ct(N+1:2*N, :); % v data

b_old = b;
b_oldx = v - s*b;

u0 = [b; v; s; theta];

% Kinetics (uses cont_df.m)

F = @(u) cont_df_h15(u, N, dx, b_old, b_oldx, M, D);

options = optimset('Jacobian', 'on', 'Display', 'iter', 'TolFun', 1e-8, ...
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
    s = u_new(2*N+1)
    theta = u_new(2*N+2)

%     plot(x, b, 'b') % space x vs biomass b
%     % axis([0 L 0 1])
%     xlabel('Space $$x$$', 'FontWeight', 'bold', 'FontSize', 24)
%     ylabel('Biomass $$b$$', 'FontWeight', 'bold', 'FontSize', 24)
%     title(['Space $$x$$ vs. Biomass $$b$$ with $$\theta$$ = ' num2str(theta) ...
%         ', $$s$$ = ' num2str(s)], 'FontWeight', 'bold', 'FontSize', 20)
               
    otherwise
        
    warning('Unexpected mode: Type 1 for debug, 2 for fsolve');
        
end