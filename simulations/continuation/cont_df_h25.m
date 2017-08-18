% To be used with continuation.m, provides function and Jacobian matrix for
% the model system detailed below with SECANT CONTINUATION: 

% b_x = v - s*b
% v_x = -(theta - v)*b^2 + b

% Zachary Singer, University of Minnesota Twin Cities, 8/10/17

%%%%%%%%%%%%%%%%%%% Initializing Function: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for fsolve set up as [F, dz]

function [f, J] = cont_dz_h25(u, N, dx, b_old, b_oldx, M, D) %[f, J]

%%%%%%%% INPUTS: %%%%%%%%%

% u     - 2*N+2 by 1 vector (z = [b; v; s; theta])
% N     - number of iteration points (scalar)
% dx    - step size (scalar)
% M     - averaging matrix (N-1 by N matrix)
% D     - first order upwind differentiation matrix (N-1 by N matrix)

%%%%%%% OUTPUTS (designed for use with fsolve): %%%%%%%%

% f     - 2*N + 2 equations, split into 2N - 2 ; 2 ; 1 ; 1 components
% J     - 2*N+2 by 2*N+2 jacobian matrix, structure shown further below

% verify dimension of u
[s1, s2] = size(u); % should be [2*N + 2, 1]
if (s1 ~= 2*N + 2 || s2 ~= 1)
    error('u needs to be size 1 by 2*N + 2');
end

% set up components
b = u(1: N);
v = u(N+1: 2*N);
s = u(2*N+1);
theta = u(2*N+2);

%%%%%%%%%%%%%%%% FIRST (N-1) + (N-1) COMPONENTS FOR b, v %%%%%%%%%%%%%%%%%%

f(1: N-1) = D*b - M*(v-s*b);
f(N: 2*N-2) = D*v - M*((v-theta).*b.^2 + b);

%%%%%%%%%%%%%%%%%%% ADDITIONAL BOUNDARY VALUES (4) %%%%%%%%%%%%%%%%%%%%%%%%

% stable eigenvalue for linearized system
lamb_pos = -0.5*s + sqrt(0.25*s^2 + 1);
lamb_pos_ds = -0.5 + 0.25*s / sqrt(0.25*s^2 + 1);

f(2*N-1) = b(1) - lamb_pos*v(1);

% v = sb, -(theta-sb)*b + 1 = 0
b_plus = theta/(2*s) + sqrt(theta^2/(4*s^2) - 1/s);
b_plus_ds = -theta/(2*s^2) + (-theta^2/(4*s^3) + 1/(2*s^2)) / ...
             (sqrt(theta^2/(4*s^2) - 1/s)); % derivative with s
b_plus_dt = 1/(2*s) + (theta/(4*s^2))/sqrt(theta^2/(4*s^2) - 1/s); % deriv w/ theta
         
% unstable left eigenvector
unstab_eval = [1 0.5*s + sqrt(0.25*s^2 + 1)];
unstab_eval_ds = [0 0.5 + 0.25*s / sqrt(0.25*s^2 + 1)]; % derivative with s

% want first component and special b_plus
init_comp = [b(N) v(N)] - [b_plus s*b_plus];
init_comp_ds = -1*[b_plus_ds b_plus + s*b_plus_ds]; % derivative with s

f(2*N) = dot(init_comp, unstab_eval);

% phase condition
f(2*N+1) = (b - b_old)' * b_oldx; % Hmm but above b_old = b...

% boundary condition
f(2*N+2) = 4*s - theta^2; 

%%%%%%%%%%%%%%%%%%%%%% JACOBIAN COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Jacobian Structure:
%
%           (b)    (v)       (s)      (theta)
%       (   -------------                    )
%       (   Jbb    Jbv   |   Jbs    |   Jbt  )
%  J =  (   Jvb    Jvv   |   Jvs    |   Jvt  )
%       (   -------------------------------- )
%       (      Jbc1      |   Jbc1s  |  Jbc1t )
%       (      Jbc2      |   Jbc2s  |  Jbc2t )
%       (               Jbc3                 )
%       (               Jbc4                 )

Jbb = D + s*M;
Jbv = -M;
Jvb = -M * spdiags(2*(v-theta).*b + ones(N,1), 0, N, N);
Jvv = D - M*spdiags(b.^2, 0, N, N);

Jbc1 = zeros(1, 2*N);
Jbc1(1, 1) = 1;             % vs. (1, N)
Jbc1(1, N+1) = -lamb_pos;   % vs. (1, 2*N)

Jbc2 = zeros(1, 2*N);
Jbc2(1, N) = unstab_eval(1);
Jbc2(1, 2*N) = unstab_eval(2);

Jbs = M*b;
Jvs = sparse(N-1, 1);

Jbc1s = lamb_pos_ds*v(1); % v(1)
Jbc2s = dot(init_comp_ds, unstab_eval) + dot(init_comp, unstab_eval_ds);

Jbt = zeros(N-1, 1);
Jvt = M*b.^2;

Jbc1t = 0;
Jbc2t = dot(-[b_plus_dt s*b_plus_dt], unstab_eval); % -

Jbc3 = zeros(1, 2*N + 2);
Jbc3(1:N) = b_oldx'; 

Jbc4 = zeros(1, 2*N + 2);
Jbc4(2*N+1) = 4;
Jbc4(2*N+2) = -2*theta;

J = [Jbb Jbv Jbs Jbt; 
     Jvb Jvv Jvs Jvt; 
     Jbc1 Jbc1s Jbc1t; 
     Jbc2 Jbc2s Jbc2t; 
     Jbc3; 
     Jbc4];
 
f = f'; % don't like row vectors

end