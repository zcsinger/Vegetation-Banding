function Df = Df(t, U, Dw, Dw_a, Db, N, dx, x, m)

w = U(1: N);
b = U(N+1: end);

g = @(w, b) w.*b.^2;

gp_b = @(w, b) 2*w.*b;
gp_w = @(w, b) b.^2;

Df = [Dw + Dw_a - spdiags(gp_w(w,b), 0, N, N)  spdiags(gp_b(w,b), 0, N, N) + m*speye(N);   
      spdiags(gp_w(w, b), 0, N, N)  Db + spdiags(gp_b(w, b), 0, N, N) - m*speye(N)];

end
