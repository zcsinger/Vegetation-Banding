function theta = maxwell(m, c, s)

% finding Maxwell condition (i.e. theta) of 
% -mb + theta*b^2 - s/(c+s) b^3 

x0 = 4*m;
a = s/(c+s);
g = @(t) @(b) -m*b + t*b.^2 - a*b.^3; % curried for parameter t and var b

rt_p = @(t) (t + sqrt(t^2 - 4*a*m))/(2*a); % upper root for g (specific)

top = @(t) integral(g(t), 0, rt_p(t)); % integral of first to last root = 0

theta = fsolve(top, x0); % find the theta that satisfies above cond.

end
