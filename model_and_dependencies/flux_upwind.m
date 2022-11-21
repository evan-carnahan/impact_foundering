function [A] = flux_upwind(q,Grid)
% author: Evan Carnahan
% date: 11/24/2019
% Description:
% This function computes the upwind flux matrix from the flux vector.
%
% Input:
% q = Nf by 1 flux vector from the flow problem.
% Grid = structure containing all pertinent information about the grid.
%
% Output:
% A = Nf by Nf matrix contining the upwinded fluxes
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> q = ones(Grid.Nf,1);S
% >> [A] = flux_upwind(q,Grid);

Nx = Grid.Nx; Ny = Grid.Ny; Nz = Grid.Nz; N = Grid.N;
Nfx = Grid.Nfx; % # of x faces
Nfy = Grid.Nfy; % # of y faces
%% One dimensinal
if ((Nx>1) && (Ny==1)) || ((Nx==1) && (Ny>1)) % 1D
   qn = min(q(1:Nfx-1),0);
   qp = max(q(2:Nfx),0);
   A = spdiags(qn,0,Nfx,Nfx-1) + spdiags(qp,-1,Nfx,Nfx-1);
%% Two dimensional
elseif (Nx>1) && (Ny>1)
    Qdp = spdiags(max(q,0),0,Grid.Nf,Grid.Nf);
    Qdn = spdiags(min(q,0),0,Grid.Nf,Grid.Nf);
%   build x-advective
    Iy = speye(Grid.Ny);
    Axp1 = spdiags(ones(Grid.Nx,1),-1,Grid.Nx+1,Grid.Nx);
    Axn1 = spdiags(ones(Grid.Nx,1),0,Grid.Nx+1,Grid.Nx);
    Axp = kron(Axp1,Iy);
    Axn = kron(Axn1,Iy);

%   build y-advictive
    Ix = speye(Grid.Nx);
    Ayp1 = spdiags(ones(Grid.Ny,1),-1,Grid.Ny+1,Grid.Ny);
    Ayn1 = spdiags(ones(Grid.Ny,1),0,Grid.Ny+1,Grid.Ny);
    Ayp = kron(Ix,Ayp1);
    Ayn = kron(Ix,Ayn1);
    %form arrays
    Ap = [Axp;Ayp];
    An = [Axn;Ayn];
    A = Qdp*Ap+Qdn*An;
   end
