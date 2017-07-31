function f = f(t, U, Dw, Dw_a, Db, N, dx, x, m)

g = @(w, b) w.*b.^2;
%
w = U(1 : N);
b = U(N + 1 : end);
%
%
f = [Dw*w + Dw_a*w - g(w, b) + m*b ; Db*b + g(w, b) - m*b];
%      
    
end