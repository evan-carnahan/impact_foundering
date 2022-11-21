function [u] = solve_lbvp(L,f,B,g,N)
% author: Evan Carnahan
% Description
% Computes the solution $u$ to the linear differential problem given by
%
% $$\mathcal{L}(u)=f \quad x\in \Omega $$
%
% with boundary conditions
%
% $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$.
%
% Input:
% L = matrix representing the discretized linear operator of size N by N,
% where N is the number of degrees of fredom
% f = column vector representing the discretized r.h.s. and contributions
% due non-homogeneous Neumann BC?s of size N by 1
% B = matrix representing the constraints arising from Dirichlet BC?s of
% sizeNcbyN
% g = column vector representing the non-homogeneous Dirichlet BC?s of size % Ncby1.
% N = matrix representing a orthonormal basis for the null-space of B and
% of size N by (N-Nc).
% Output:
% u = column vector of the solution of size N by 1
if isempty(B) % no constraints
    u = L\f;
else
    u_p = spalloc(length(f),1,length(g));
    u_p = B'*((B*B')\g);
    u_0 = N*((N'*L*N)\(N'*(f-(L*u_p))));
    u = u_0+u_p;
end
