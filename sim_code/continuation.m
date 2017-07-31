% Solve a particular vegetative system through continuation methods

% Zachary Singer, University of Minnesota Twin Cities, 7/28/17

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scaled system: 
% b_x = v - s*b
% v_x = -(theta - v)*b^2 + b

% theta - varied parameter (want to find theta in terms of s)

%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
% set(gcf, 'Position', [100, 800, 800, 800])

m = input('Enter a mode (1 for normal execution, 2 for Jacobian debug): ');

% Step Size Constants
L = 25;     % conservative placeholder
dx = 0.2;
x = (-L: dx: L)'; % spacing of dx up to L, transpose
N_tot = size(x');
N = N_tot(2); % how many time steps

% Set up matrix form
e = ones(N, 1);
D = spdiags([-e e], 0: 1, N-1, N) / dx;

M = spdiags([e e], 0: 1, N-1, N) / 2;

% Parameters in equation
% theta = .7;
% s = -1;
theta = 3;
s = 1;

% Placeholder Phase Condition
cnst = 5;
b_amp = theta/(2*s) + sqrt(theta^2/(4*s^2) - 1/s);
b_old = b_amp ./ (1+exp(cnst*x));
b_oldx =  -b_amp*cnst * exp(cnst*x) ./ (exp(cnst*x) + 1).^2; 
b_oldxx = b_amp*cnst^2 * exp(cnst*x) .* (exp(cnst*x) - 1) ...
            ./ (exp(cnst*x) + 1).^3;

u0 = [b_old; b_oldx; s];

F = @(u) cont_df(u, N, dx, theta, b_old, b_oldx, b_oldxx, M, D);

options = optimset('Jacobian', 'on','Display','iter','TolFun',1e-6, ...
          'TolX',1e-6,'MaxIter',30,'Algorithm','trust-region-reflective');

%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING/TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch m % mode
    
    case 2

    ord = input('Enter an order of magnitude (negative): ');
    disp(['Testing Jacobian with order 10^' num2str(ord) '...'])
    du = 10^ord * rand(2*N + 1, 1); % random points, order 'ord'
    [F0, J0] = F(u0);
    [F1, J1] = F(u0 + du);
    err_u = (F1 - F0 - J0*du);
    disp(['Error: ' num2str(norm(err_u))]);
    plot(err_u)


%%%%%%%%%%%%%%%%%%% Initial use of fsolve: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 1

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

    otherwise
        
    warning('Unexpected mode: Type 1 for normal execution, 2 for debug.');
        
end
