% Solves nonlinear system to find new dz direction

% Zachary Singer, University of Minnesota - Twin Cities 8/1/17

function f = dz_solver(dz, dz_old, Jz, N)

% Inputs:
% dz     - direction (1 by 2*N+2 vector) 
% dz_old - old direction (1 by 2*N+2 vector)
% Jz     - Jacobian (2*N+1 by 2*N+2 matrix)

f(1:2*N+1) = Jz * dz';

f(2*N+2) = dot(dz, dz_old) - 1;

end