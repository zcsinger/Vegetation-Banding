function F = veg_fb(z, b_old, b_oldx, D, D2, N, m, c, tau, theta, dx)

% Model: b_xx + s*b_x - m*b + theta*b^2 - s/(c+s)*b^3 - tau/(c+s)*b_x*b^2

thres = 5; % threshold around center
b = z(1:N);
s = z(N+1);

F(1:N) = D2*b + s*D*b - m*b + theta*b.^2 - s/(c+s)*b.^3 - tau/(c+s)*(D*b).*b.^2;

int = floor(N/2) - thres : floor(N/2) + thres;
A = z(int)' - b_old(dx*int);
B = b_oldx(dx*int);

F(N+1) = dot(A, B);

end