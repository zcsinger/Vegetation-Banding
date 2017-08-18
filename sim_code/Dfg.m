function Df=Df(t,U,Dw,Dw_a,Db,N,dx,x,gam,a)

w=U(1:N);b=U(N+1:end);

g=@(b) -b.*(1-b).*(b-a);
gp=@(b) -( (1-b).*(b-a)+b.*(1-b)-b.*(b-a) );

Df=[Dw + Dw_a - gam*speye(N)    spdiags(gp(b),0,N,N);   
    +gam*speye(N)  Db-spdiags(gp(b),0,N,N)];

end
