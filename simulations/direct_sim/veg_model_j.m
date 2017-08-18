% finds jacobian evaluation for model below, used in veg_model_80317.m: 

% Zachary Singer, University of Minnesota - Twin Cities, 8/3/17

% b_t = b_xx + w*b^2 - b 
% w_t =      - w*b^2 + b + cw_x

function J = veg_model_j(t, V, D1, D2, N, c)

b = V(1: N);
w = V(N+1: 2*N);

% F(1: N) = D2*b + w .* b.^2 - b;
% F(N+1 : 2*N) = -w .* b.^2 + b + c*D*w;

J_bb = D2 + spdiags(2*w .* b - ones(N, 1), 0, N, N);
J_bw = spdiags(b.^2, 0, N, N);
J_wb = spdiags(-2*w .* b + ones(N, 1), 0, N, N);
J_ww = spdiags(-b.^2, 0, N, N) + c*D1;

J = [J_bb J_bw; J_wb J_ww];

end