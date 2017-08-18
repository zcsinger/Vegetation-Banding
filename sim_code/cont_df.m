function [f, J] = cont_df(u, N, dx, theta, b_old, b_oldx, b_oldxx, M, D)

% f : 2*N + 1 equations, split into 2N - 2 ; 2 ; 1 
% b_x = v_x - s*b
% v_x = (v - theta)*b^2 + b

% u = [b; v; s] is vector of size 2*N + 1

% verify dimension of u
[s1, s2] = size(u); % should be [1, 2*N + 1]
if (s1 ~= 2*N + 1 || s2 ~= 1)
    error('u needs to be size 1 by 2*N + 1');
end

% set up components
b = u(1: N);
v = u(N+1: 2*N);
s = u(2*N+1);

%%%%%%%%%% FIRST (N-1) + (N-1) COMPONENTS %%%%%%%%%%%

f(1: N-1) = D*b - M*(v-s*b);
f(N: 2*N-2) = D*v - M*((v-theta).*b.^2 + b);

% vectorized expressions, potential sign error
% j = 1: N-1;
% f(j) = (b(j+1) - b(j)) / dx - (v(j) - s*b(j)) / 2 - (v(j+1) - s*b(j+1)) / 2;
% 
% k = N: 2*N-2; % want to set f(N : 2*N-2)
% f(k) = (v(j+1) - v(j)) / dx - (v(j+1) .* b(j+1).^2 - theta*b(j+1).^2 + ...
%         b(j+1)) / 2 - (v(j) .* b(j).^2 - theta*b(j).^2 + b(j)) / 2;

%%%%%%%%%% ADDITIONAL BOUNDARY VALUES (3) %%%%%%%%%%%

% stable eigenvalue for linearized system
lamb_neg = -0.5*s - sqrt(0.25*s^2 + 1);
lamb_neg_ds = -0.5 - 0.25*s / sqrt(0.25*s^2 + 1);

f(2*N-1) = b(N) - lamb_neg*v(N);

% v = sb, -(theta-sb)*b + 1 = 0
b_plus = theta/(2*s) + sqrt(theta^2/(4*s^2) - 1/s);
b_plus_ds = -theta/(2*s^2) + (-theta^2/(4*s^3) + 1/(2*s^2)) / ...
             (sqrt(theta^2/4*s^2 - 1/s)); % derivative with s

% stable left eigenvector
stab_eval = [1 0.5*s - sqrt(0.25*s^2 + 1)];
stab_eval_ds = [0 0.5 - 0.25*s / sqrt(0.25*s^2 + 1)]; % derivative with s

% stab_eval = [1 0];
init_comp = [b(1) v(1)] - [b_plus s*b_plus];
init_comp_ds = -1*[b_plus_ds b_plus + s*b_plus_ds]; % derivative with s

f(2*N) = dot(init_comp, stab_eval);

% phase condition
f(2*N+1) = (b - b_old)'*b_oldx;

%%%%%%%%%%%%%%%%%%%%%% JACOBIAN COMPONENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jbb = D + s*M;
Jbv = -M;
Jvb = -M * spdiags(2*(v-theta).*b + ones(N,1), 0, N, N);
Jvv = D - M*spdiags(b.^2, 0, N, N);

Jbc1 = zeros(1, 2*N);
Jbc1(1, N) = 1;
Jbc1(1, 2*N) = -lamb_neg;

Jbc2 = zeros(1, 2*N);
Jbc2(1, 1) = stab_eval(1);
Jbc2(1, N+1) = stab_eval(2);

Jbs = M*b;
Jvs = sparse(N-1, 1);

Jbc1s = lamb_neg_ds*v(N); 
Jbc2s = dot(init_comp_ds, stab_eval) + dot(init_comp, stab_eval_ds);

Jbc3 = zeros(1, 2*N + 1);
Jbc3(1:N) = b_oldx'; % *(dot(ones(N, 1) - b_oldx, b_oldx) + dot(b-b_oldx, b_oldxx));


J = [Jbb Jbv Jbs; Jvb Jvv Jvs; Jbc1 Jbc1s; Jbc2 Jbc2s; Jbc3];
f = f'; % i don't row vector

%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING/TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jacobian Sizes (debugging)
% disp(['size of Jbb: ' num2str(size(Jbb))])
% disp(['size of Jbv: ' num2str(size(Jbv))])
% disp(['size of Jbs: ' num2str(size(Jbs))])
% disp(['size of Jvb: ' num2str(size(Jvb))])
% disp(['size of Jvv: ' num2str(size(Jvv))])
% disp(['size of Jvs: ' num2str(size(Jvs))])
% disp(['size of Jbc1: ' num2str(size(Jbc1))])
% disp(['size of Jbc1s: ' num2str(size(Jbc1s))])
% disp(['size of Jbc2: ' num2str(size(Jbc2))])
% disp(['size of Jbc2s: ' num2str(size(Jbc2s))])
% disp(['size of Jbc3: ' num2str(size(Jbc3))])

end
