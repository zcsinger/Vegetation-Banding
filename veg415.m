close all;

%%%%%%%%%%%%%
% integrator for coupled transport equations with turning rates in fg.m, derivatives of turning rates in Dfg.m

% w_t = dw*w_xx + (c*w_x) + g(b) - gam*w
% b_t = db*b_xx - g(b) + gam*w

% also option of depositing mass w, b into the system; see fg.m, "trigger" (in this case change initial condition to U=0*U

% spatial discritization via upwind derivatives, 1st, 2nd, or 3rd order

% time stepping with explicit Euler, step size in terms of fraction of dx (Dfg not needed)

% alternate time stepping with ode15s


% also plots ODE bifurcation diagram when tr(u,v)=u(1+v^2)/(1+gamma v^2)
% needs gamma as parameters



%%%% initialize grid %%%%%%%%%%%%%%%%%%%

L = 400;
dx = .3;
x = (dx:dx:L)';
N = size(x);
N = N(1);
dw = 1;
c = 0.5;
db = .1;
a = .5;
%%%%  upwind derivatives %%%%%%%%%%%%%%%%%%%
e = ones(N,1);

%%% 1st order approximation, introducing viscosity rather than dispersion
%  Dp=spdiags([-e e],-1:0,N,N);Dp(1,N)=-1; sf=0.8;

% diffusion of water
Dw = dw*spdiags([e -2*e e],-1:1,N,N)/dx^2;
% Neumann Boundary Conditions
Dw(1,1) = -dw/dx^2;
Dw(N,N) = -dw/dx^2;

% advection of water
Dw_a = c*spdiags([-e e], 0:1, N, N) / dx; % check this
% upwind is [e -e] vs. downwind [-e e]
% diags could be either 0:1 or -1:0
% Dw_a(1, 1) = -1/dx;
Dw_a(N, N) = 0;

% diffusion of biomass
Db = db*spdiags([e -2*e e],-1:1,N,N)/dx^2;
% Neumann Boundary Conditions
Db(1,1) = -db/dx^2;
Db(N,N) = -db/dx^2;

%%%% nonlinearity, kinetics %%%%%%%%%%%%%%%%%%%%%%%%%%

gam = 0.1;
f0 = @(t,V) fg(t,V,Dw,Dw_a,Db,N,dx,x,gam,a);
Df0 = @(t,V) Dfg(t,V,Dw,Dw_a,Db,N,dx,x,gam,a);

%%%%  time stepping parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('RelTol',1e-3,'AbsTol',1e-3,'Jacobian',Df0);

tmax = 2*L;
dt = 2;
tspan =[0 dt];

%%%%  initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

al = a; % u component  of initial data
% symmetric initial data (comment one out)
  w = 0*ones(N,1);
  % w(50:100) = .3; test for 
  b = 1*ones(N,1);  
w(1:10) = w(1:10)+.01; % perturbation
b(1 :floor(N/4)) = 0;  % creating interface in center
b(floor(3*N/4) : N) = 0;


% concatenate initial data to form long vector
U = [w ; b];
%  U=0*U; %  uncomment when turning on trigger, choose trigger amplitude amp to a then in fg.m

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
for j = 1:tmax/dt
    [T,UALL] = ode15s(f0,tspan+j*dt,U,options);
    U=UALL(end,:)';
   
    w = U(1:N);
    b = U(N+1:2*N);
    % PLOTTING ROUTINE %%%%
    if plottrue
      % set(0, 'CurrentFigure', h);
      subplot(3,1,1)
      plot(x,w,'r')    
        axis([0 L -.40  a])
        xlabel('Space x')
        ylabel('Water w')
        title(sprintf('Space vs. Water with c = %.2f', c))
      subplot(3,1,2)
      plot(x,b,'b')
        axis([0 L -.40  2*a])
        xlabel('Space x')
        ylabel('Biomass b')
        title(['Space vs. Biomass with c = ' num2str(c) ' and time ' num2str(j*dt)])
      subplot(3,1,3)
      plot(b,w,'g')
        axis([0 1 -.5 .5])
        xlabel('Biomass b')
        ylabel('Water w')
        title(sprintf('Biomass vs. Water with c = %.2f', c))
      drawnow
%        pause
    end
   
    % END PLOTTING ROUTINE %%%%
end


