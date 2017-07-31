function F = veg_fb(z, b_old, b_oldx, D, D2, N, m, c, tau, theta, dx)

% Skeletal Model: b_xx + sb_x + b - b^3 = 0
% Explicit Solution : tanh(x/sqrt(2))

thres = 5; % threshold around center
b = z(1:N);
s = z(N+1);

F(1:N) = D2*b + b - b.^3 + s*D*b; % turn on/off speed

int = floor(N/2) - thres : floor(N/2) + thres;
A = z(int)' - b_old(dx*int);
B = b_oldx(dx*int);

F(N+1) = dot(A, B);

end