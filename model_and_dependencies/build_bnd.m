function [B,N,fn] = build_bnd(Param,Grid,I)
   % author: Evan Carnahan
   % date: 9/20/19
   % Description:
   % This function computes the operators and r.h.s vectors for both Dirichlet
   % and Neumann boundary conditions.
   %
   % Input:
   % Grid = structure containing all pertinent information about the grid.
   % Param = structure containing all information about the physical problem
   %         in particular this function needs the fields
   %         Param.dof_dir = Nc by 1 column vector containing
   %                         the dof?s of the Dirichlet boundary.
   %         Param.dof_neu = column vector containing
   %                         the dof?s of the Neumann boundary.
   %         Param.qb      = column vector of prescribed fluxes on Neuman bnd.
   % I = identity matrix in the full space
   %% If Neumann condtions are not set explicitly, make them natural
if ~isfield(Param,'dof_neu');   Param.dof_neu   = []; end
if ~isfield(Param,'dof_f_neu'); Param.dof_f_neu = []; end
if ~isfield(Param,'qb');        Param.qb        = []; end
   
   B = I([Param.dof_dir],:);
   N = I; 
   N(:,[Param.dof_dir]) = [];
   
   %Added 9/30
%    fn = spalloc(Grid.Nx,1,length(Param.qb));
%    fn(Param.dof_neu) = Param.qb.*Grid.A(Param.dof_neu)./Grid.V(Param.dof_neu);
%    fn = [];
if isempty(Param.dof_neu)
    fn = spalloc(Grid.N,1,0);
else
    fn = spalloc(Grid.N,1,length(Param.dof_neu));
    fn(Param.dof_neu) = Param.qb.*Grid.A(Param.dof_f_neu)./Grid.V(Param.dof_neu);
end
   % Output:
   % B = Nc by N matrix of the Dirichlet constraints
   % N = (N-Nc) by (N-Nc) matrix of the nullspace of B
   % fn = N by 1 r.h.s. vector of Neuman contributions
   %
   % Example call:
   % >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
   % >> Grid = build_grid(Grid);
   % >> [D,G,I]=build_ops(Grid);
   % >> Param.dof_dir   = Grid.dof_xmin;    % identify cells on Dirichlet bnd
   % >> Param.dof_f_dir = Grid.dof_f_xmin;  % identify faces on Dirichlet bnd
   % >> Param.dof_neu   = Grid.dof_xmax;    % identify cells on Neumann bnd
   % >> Param.dof_f_neu = Grid.dof_f_xmax;  % identify cells on Neumann bnd
   % >> Param.qb = 1;                   % set bnd flux
   % >> [B,N,fn] = build_bnd(Param,Grid,I);