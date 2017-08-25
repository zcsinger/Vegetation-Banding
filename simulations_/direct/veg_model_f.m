% finds function evaluation for model below, used in veg_model_80317.m: 

% Zachary Singer, University of Minnesota - Twin Cities, 8/3/17

% b_t = b_xx + w*b^2 - b 
% w_t =      - w*b^2 + b + cw_x

% (matrix form)
% 0 = D2*b + w .* b.^2 - b;
% 0 =      - w .* b.^2 + b + c*D1*w;

function F = veg_model_f(t, V, D1, D2, N, dx, c, w_plus, b_minus)

w_minus = 1 / b_minus;

b = V(1: N);
w = V(N+1: 2*N);

F(1: N) = D2*b + w .* b.^2 - b;

F(1) = F(1) + b_minus / dx^2;
% F(N) = F(N) + b_plus / dx^2; % biomass boundary condition on right

F(N+1 : 2*N) = -w .* b.^2 + b + c*D1*w;

F(N+1) = F(N+1) + c * w_minus / dx;
F(2*N) = F(2*N) + c * w_plus / dx; % water boundary condition

F = F'; % column vector

end