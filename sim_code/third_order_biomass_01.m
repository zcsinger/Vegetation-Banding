% solve third order ODE below:
% db*dw*b_xxx - db*c*b_xx - (dw*gp+db*gam)*b_x+c*g+gam*psi

% constants
L = 100;
dx = .2;
x = (dx:dx:L)';
N = size(x);
N = N(1);
c = 0.1;
a = 0.5;
gam = 0.5;
psi =  1; % technically a variable
dw = 1;
db = 0.01; % db << dw


%%%%  time stepping parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = @(t,V) fg2(t, V, a, gam, psi, c, db, dw);
% Df0 = @(t,V) Dfg(t,V,Dw,Dw_a,Db,N,dx,x,gam,a);

% options = odeset('RelTol',1e-3,'AbsTol',1e-3,'Jacobian',Df0);

tmax = 2*L;
dt = 2;
t_end = 10 ;
tspan =[0 t_end];

% initial condition
x0 = [1 0 0];

[t, b] = ode15s(f0, tspan, x0);
figure
subplot(2,1,1)
plot(b(:, 1), b(:, 2))
xlabel('b')
ylabel('b_x')
title(sprintf('phase portrait for third order b_{xxx} equation, psi = %.2f', psi))
subplot(2,1,2)
plot(t, b(:, 1))
xlabel('Space x')
ylabel('Biomass b')
title(sprintf('time versus b for third order b_{xxx} equation, psi = %.2f', psi))
% subplot(3,1,3)
% plot3(b(:, 1), b(:, 2), b(:, 3))
% xlabel('Biomass b')
% ylabel('b_x')
% zlabel('b_{xx}')
% title(sprintf('3D phase space for third order b_{xxx} equation, psi = %.2f', psi))

%%%%  initial condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% al = a; % u component  of initial data
% symmetric initial data (comment one out)
  % w = 0*ones(N,1);
  % b = al*ones(N,1);
% w(1:10) = w(1:10)+.01;

% concatenate initial data to form long vector
% U = [w;b];
% U=0*U; %  uncomment when turning on trigger, choose trigger amplitude amp to a then in fg.m


% % plotting initialize
% plottrue=1;
% savetrue=1;
% % initialize plotting
% if plottrue 
%   h = figure(1);
%     drawnow
% end


%%%% main time stepping routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j = 1 : tmax/dt
%     [T,b_all] = ode15s(f0,tspan+j*dt,b); % options   
%     b = b_all(end,:)';
%     % PLOTTING ROUTINE %%%%
%     if plottrue
%         plot(x, b);
%         hold on
% %       pause
%     end
%    
%     % END PLOTTING ROUTINE %%%%
% end


