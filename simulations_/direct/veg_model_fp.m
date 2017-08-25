% finds function evaluation for model below, used in veg_model_hom.m: 
% periodic boundary conditions

% Zachary Singer, University of Minnesota - Twin Cities, 8/16/17

% b_t = b_xx + w*b^2 - b 
% w_t =      - w*b^2 + b + cw_x

% (matrix form)
% 0 = D2*b + w .* b.^2 - b;
% 0 =      - w .* b.^2 + b + c*D1*w;

% Difference : Periodic Boundary Conditions

function F = veg_model_fp(t, V, D1, D2, N, dx, c, w_plus)

b = V(1: N);
w = V(N+1: 2*N);

F(1: N) = D2*b + w .* b.^2 - b;
F(N+1 : 2*N) = -w .* b.^2 + b + c*D1*w;
F(N) = F(1);        % periodic bc in biomass
F(2*N) = F(N+1);    % periodic bc in water
F = F'; % column vector

end