function [Kd] = comp_mean_corners(K,p,Grid)
% author: Evan Carnahan
% Description:
% Takes coefficient field, K, defined at the cell centers and computes the
% mean specified by the power, p and returns vector of viscosities, Kd.
%
% Input:
% K = Ny by Nx matrix of cell centered values
% p = power of the generalized mean
%       1 (arithmetic mean)
%      -1 (harmonic mean)
% Grid = structure containing information about the grid.

mean_k = nan(Grid.Ny+1,Grid.Nx+1);

% use 4-point stencil to calculate internal corner property
mean_k(2:Grid.Ny,2:Grid.Nx) = ((K(1:Grid.Ny-1,1:Grid.Nx-1).^p+K(2:Grid.Ny,2:Grid.Nx).^p...
    +K(1:Grid.Ny-1,2:Grid.Nx).^p+K(2:Grid.Ny,1:Grid.Nx-1).^p)/4).^(1/p);

%Find edge corners with 2-point stencil
mean_k(Grid.Ny+1,2:Grid.Nx) = ((K(Grid.Ny,1:Grid.Nx-1).^p+K(Grid.Ny,2:Grid.Nx).^p)/2).^(1/p);
mean_k(1,2:Grid.Nx) = ((K(1,2:Grid.Nx).^p+K(1,1:Grid.Nx-1).^p)/2).^(1/p);
mean_k(2:Grid.Ny,Grid.Nx+1) = ((K(2:Grid.Ny,Grid.Nx).^p+K(1:Grid.Ny-1,Grid.Nx).^p)/2).^(1/p);
mean_k(2:Grid.Ny,1) = ((K(1:Grid.Ny-1,1).^p+K(2:Grid.Ny,1).^p)/2).^(1/p);

%return vector of computed quantity at corners
KdVec = mean_k(:);
Kd = KdVec;

   