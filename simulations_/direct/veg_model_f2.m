% finds function evaluation for model below, used in veg_model_80817.m: 

% Zachary Singer, University of Minnesota - Twin Cities, 8/8/17

% b_t = b_xx + w*b^2 - b 
% w_t =      - w*b^2 + b + cw_x

% (matrix form)
% 0 = D2*b + w .* b.^2 - b;
% 0 =      - w .* b.^2 + b + c*D1*w;

% Difference : Boundary condition on left

function F = veg_model_f2(t, V, D1, D2, N, dx, c, b_plus)

w_plus = 1 / b_plus;

b = V(1: N);
w = V(N+1: 2*N);

F(1: N) = D2*b + w .* b.^2 - b;
F(N+1 : 2*N) = -w .* b.^2 + b + c*D1*w;
F(N) = F(N) + b_plus / dx^2; % biomass boundary condition on right
F(2*N) = F(2*N) + c * w_plus / dx; % water boundary condition
F = F'; % column vector

end