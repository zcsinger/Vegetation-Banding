function f=f(t,U,Dw,Dw_a,Db,N,dx,x,gam,a)

g = @(b) -b.*(1-b).*(b-a);
%
w = U(1:N);b=U(N+1:2*N);
%
%
f = [Dw*w + Dw_a*w + g(b) - gam*w; Db*b - g(b) + gam*w];
%    
%    
%  % trigger
%      amp=5;del=0.3;s=0.2;
%      source=amp*s/del/pi*(sech((x-s*t)/del)+sech((x+s*t-max(x))/del));
%      f=f+[1*source;1*source];
%      
    
end
  