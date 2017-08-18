% To be used with continuation.m, provides function and Jacobian matrix for
% the model system detailed below with SECANT CONTINUATION: 

% b_x = v - s*b
% v_x = -(theta - v)*b^2 + b

% Zachary Singer, University of Minnesota Twin Cities, 8/1/17

%%%%%%%%%%%%%%%%%%% Initializing Function: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for fsolve set up as [F, dz]

function [f, J, Dz] = cont_dz_h2(z, z_old, N, dx, ds, dz, b_old, b_oldx, M, D)

%%%%%%%% INPUTS: %%%%%%%%%

% z     - 2*N+2 by 1 vector (z = [b; v; s; theta])
% z_old - old data (also 2*N+2 by 1 vector)
% N     - number of iteration points (scalar)
% dx    - step size (scalar)
% ds    - arclength step size (scalar ~ 10^-2)
% dz    - direction vector (2*N+2 by 1 vector)
% b_old - old b vector (N by 1 vector)
% b_oldx- old b_x/v vector (N by 1 vector)
% M     - averaging matrix (N-1 by N matrix)
% D     - first order upwind differentiation matrix (N-1 by N matrix)

%%%%%%% OUTPUTS (designed for use with fsolve): %%%%%%%%

% f     - 2*N + 2 equations, split into 2N - 2 ; 2 ; 1 ; 1 components
% J     - 2*N+2 by 2*N+2 jacobian matrix, structure shown further below

% verify dimension of u
[s1, s2] = size(z); % should be [2*N + 2, 1]
if (s1 ~= 2*N + 2 || s2 ~= 1)
    error('z needs to be size 1 by 2*N + 2');
end

% set up components
b = z(1: N);
v = z(N+1: 2*N);
s = z(2*N+1);
theta = z(2*N+2);

%%%%%%%%%%%%%%%% FIRST (N-1) + (N-1) COMPONENTS FOR b, v %%%%%%%%%%%%%%%%%%

f(1: N-1) = D*b - M*(v-s*b);
f(N: 2*N-2) = D*v - M*((v-theta).*b.^2 + b);

% vectorized expressions, not as efficient but possibly easier to read
% j = 1: N-1;
% f(j) = (b(j+1) - b(j))/dx - (v(j) - s*b(j))/ 2 - (v(j+1) - s*b(j+1))/ 2;
% 
% k = N: 2*N-2; % want to set f(N : 2*N-2)
% f(k) = (v(j+1) - v(j)) / dx - (v(j+1) .* b(j+1).^2 - theta*b(j+1).^2 + ...
%         b(j+1)) / 2 - (v(j) .* b(j).^2 - theta*b(j).^2 + b(j)) / 2;

%%%%%%%%%%%%%%%%%%% ADDITIONAL BOUNDARY VALUES (4) %%%%%%%%%%%%%%%%%%%%%%%%

% stable eigenvalue for linearized system
lamb_pos = -0.5*s + sqrt(0.25*s^2 + 1);
lamb_pos_ds = -0.5 + 0.25*s / sqrt(0.25*s^2 + 1);

f(2*N-1) = b(1) - lamb_pos*v(1);

% v = sb, -(theta-sb)*b + 1 = 0
b_plus = theta/(2*s) + sqrt(theta^2/(4*s^2) - 1/s);
b_plus_ds = -theta/(2*s^2) + (-theta^2/(4*s^3) + 1/(2*s^2)) / ...
             (sqrt(theta^2/(4*s^2) - 1/s)); % derivative with s
b_plus_dt = 1/(2*s) + theta/sqrt(theta^2/(4*s^2) - 1/s); % deriv w/ theta
         
% stable left eigenvector
unstab_eval = [1 0.5*s + sqrt(0.25*s^2 + 1)];
unstab_eval_ds = [0 0.5 + 0.25*s / sqrt(0.25*s^2 + 1)]; % derivative with s

% want first component and special b_plus
init_comp = [b(N) v(N)] - [b_plus s*b_plus];
init_comp_ds = -1*[b_plus_ds b_plus + s*b_plus_ds]; % derivative with s

f(2*N) = dot(init_comp, unstab_eval);

% phase condition
f(2*N+1) = (b - b_old)' * b_oldx;

% want perpendicular to direction chosen
f(2*N+2) = dot(dz', z-(z_old + dz'));

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
%       (            Jbc3             (Jbc3t))
%       (               Jbcz                 )

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
Jbc2t = -[b_plus_dt s*b_plus_dt] * unstab_eval';

Jbc3 = zeros(1, 2*N + 2);
Jbc3(1:N) = b_oldx'; 

Jbcz = dz;

J = [Jbb Jbv Jbs Jbt; 
     Jvb Jvv Jvs Jvt; 
     Jbc1 Jbc1s Jbc1t; 
     Jbc2 Jbc2s Jbc2t; 
     Jbc3; 
     Jbcz];
 
f = f'; % don't like row vectors


%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING/TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jacobian Sizes (debugging)
% disp(['size of Jbb: ' num2str(size(Jbb))])
% disp(['size of Jbv: ' num2str(size(Jbv))])
% disp(['size of Jbs: ' num2str(size(Jbs))])
% disp(['size of Jbt: ' num2str(size(Jbt))])
% disp(['size of Jvb: ' num2str(size(Jvb))])
% disp(['size of Jvv: ' num2str(size(Jvv))])
% disp(['size of Jvs: ' num2str(size(Jvs))])
% disp(['size of Jvt: ' num2str(size(Jvt))])
% disp(['size of Jbc1: ' num2str(size(Jbc1))])
% disp(['size of Jbc1s: ' num2str(size(Jbc1s))])
% disp(['size of Jbc1t: ' num2str(size(Jbc1t))])
% disp(['size of Jbc2: ' num2str(size(Jbc2))])
% disp(['size of Jbc2s: ' num2str(size(Jbc2s))])
% disp(['size of Jbc2t: ' num2str(size(Jbc2t))])
% disp(['size of Jbc3: ' num2str(size(Jbc3))])
% disp(['size of Jbcz: ' num2str(size(Jbcz))])

end