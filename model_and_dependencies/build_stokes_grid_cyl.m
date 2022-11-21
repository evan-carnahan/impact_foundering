function [Grid] = build_stokes_grid_cyl(Gridp)
    %build_stokes_grid.m
    %Authors: Evan Carnahan, Marc Hesse
  
    %build pressure grid
    Gridp = build_grid(Gridp);
    
    %build x-vel grid
    Gridx.xmin = Gridp.xmin-Gridp.dx/2; Gridx.xmax = Gridp.xmax+Gridp.dx/2; Gridx.Nx = Gridp.Nx+1;
    Gridx.ymin = Gridp.ymin; Gridx.ymax = Gridp.ymax; Gridx.Ny = Gridp.Ny;
    Gridx.geom = Gridp.geom;
    Gridx = build_grid(Gridx);

    %build y-vel grid
    Gridy.xmin = Gridp.xmin; Gridy.xmax = Gridp.xmax; Gridy.Nx = Gridp.Nx;
    Gridy.ymin = Gridp.ymin-Gridp.dy/2; Gridy.ymax = Gridp.ymax+Gridp.dy/2; Gridy.Ny = Gridp.Ny+1;
    Gridy.geom = Gridp.geom;
    Gridy = build_grid(Gridy);

    Grid.p = Gridp;
    Grid.x = Gridx;
    Grid.y = Gridy;
    
    %Helpful summary features
    Grid.N = Grid.p.Nf+Grid.p.N;
    

