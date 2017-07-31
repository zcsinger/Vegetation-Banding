function Dv = veg_dfb(z, b_old, b_oldx, D, D2, N, m, c, tau, theta, dx)

b = z(1:N);
s = z(N+1);

% v(1:N) = D2*b + s*D*b - m*b + theta*b.^2 - s/(c+s)*b.^3 - tau/(c+s)*(D*b).*b.^2;

b2_spdiags = spdiags(b.^2, 0, N, N);
b_spdiags = spdiags(b, 0, N, N);

Dv = D2 + s*D - m*speye(N) + 2*theta* b_spdiags - ...
            3*s/(c+s)* b2_spdiags - ...
            tau/(c+s)* (D.* b2_spdiags + 2*(D*. b_spdiags).* b_spdiags);

A = b - b_old;
B = b_oldx;

% v(N+1) = dot(A, B);

end
