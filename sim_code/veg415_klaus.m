close all;
clc
clear all

%%%%%%%%%%%%%

% integrator for coupled transport equations with turning rates in fgkl.m, 
% derivatives of turning rates in Dfgkl.m

% altered Klausmeier model
% w_t = -wb^2 + mb + cw_x + dw*w_xx
% b_t = wb^2 - mb +db*b_xx

%%%% initialize grid and constants %%%%%%%%%%%%%%%%%%%

L = 400;    % grid length
dx = .3;    % grid spacing
x = (dx:dx:L)';
N = size(x);
N = N(1);
dw = 1;     % water diffusivity
c = .75;    % advection constant
db = .01;   % biomass diffusivity with db << dw
m = .5;      % plant loss/growth constant

%%%%%%%%%%%%%%%%%%  upwind derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = ones(N,1);

% diffusion of water
Dw = 0*dw*spdiags([e -2*e e],-1:1,N,N)/dx^2;
% Neumann Boundary Conditions
% Dw(1,1) = -dw/dx^2;
% Dw(N,N) = -dw/dx^2;

% advection of water
Dw_a = c*spdiags([-e e], 0:1, N, N) / dx; % upwind 
% upwind is [e -e] vs. downwind [-e e]

Dw_a(N, N) = 0;

% diffusion of biomass
Db = db*spdiags([e -2*e e],-1:1,N,N)/dx^2;
% Neumann Boundary Conditions
Db(1,1) = -db/dx^2;
Db(N,N) = -db/dx^2;

%%%% nonlinearity, kinetics %%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = @(t,V) fgkl(t,V,Dw,Dw_a,Db, N,dx, x, m);
Df0 = @(t,V) Dfgkl(t,V,Dw,Dw_a,Db, N,dx, x, m);

%%%%  time stepping parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('RelTol',1e-3,'AbsTol',1e-3,'Jacobian',Df0);

tmax = 200;
dt = 2;
tspan =[0 dt];

%%%%  initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% additional parameters for finding speed of biomass
b0 = 0.05;
w0 = -m/b0;

b_start = 3*N/4;
b_old = b_start * dx;
b_new = 0;
speed = 0;
mean_speed = 0;
speed_vec = [];
thres = 3;
old_time = 1;

%%%%%%%%%%%%%%%%%%% Setting up Interface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 0*ones(N,1);
b = 3*ones(N,1);  

% b(1 : floor(N/4)) = b0; % initial part of interfact
% w(1 : floor(N/4)) = w0;
% b(floor(b_start) : N) = b0;
% w(floor(b_start) : N) = w0;

% simple interface (blob in middle half)

b(1 : floor(N/4)) = 0;
b(floor(3*N/4) : N) = 0;
w(1 : floor(N/4)) = 3;
w(floor(3*N/4) : N) = 3;
% w = 3 - b;

% concatenate initial data to form long vector
U = [w ; b];

%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting initialize
plottrue=1;
savetrue=1;
% initialize plotting
if plottrue 
  h=figure(1);
    drawnow
end

%%%% main time stepping routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : tmax/dt
    [T, UALL] = ode15s(f0, tspan+j*dt, U, options);
    U = UALL(end,:)';
   
    w = U(1:N);
    b = U(N+1:2*N);
    
    b_new = fj(b, N, dx)
    if b_new - b_old > thres
        speed = (b_new - b_old)/(j - old_time);
        speed_vec = [speed_vec speed];
        mean_speed = mean(speed_vec);
        old_time = j;
        b_old = b_new;
    end
        
    % PLOTTING ROUTINE %%%%
    if plottrue
      % set(0, 'CurrentFigure', h);
      subplot(2,1,1)
      plot(x,w,'r')    
        %axis([0 L -.5  .5])
        xlabel('Space x')
        ylabel('Water w')
        title(['Space vs. Water with speed = ' num2str(speed)]) % placeholder
      subplot(2,1,2)
      plot(x,b,'b')
        %axis([0 L -.40  2])
        xlabel('Space x')
        ylabel('Biomass b')
        title(['Space vs. Biomass with c = ' num2str(c) ' and time ' num2str(j*dt)])
%       subplot(3,1,3)
%       plot(b,w,'g')
%         axis([0 1 -.5 .5])
%         xlabel('Biomass b')
%         ylabel('Water w')
%         title(sprintf('Biomass vs. Water with c = %.2f', c))
      drawnow
%        pause
    end
   
    % END PLOTTING ROUTINE %%%%
end


