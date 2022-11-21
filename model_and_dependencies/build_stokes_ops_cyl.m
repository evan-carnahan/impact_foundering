function [D,Edot,Dp,Gp,I,Gyy]=build_stokes_ops_cyl(Grid)
    % build_stokes_ops.m
    % authors: Evan Carnahan, Marc Hesse
    [Dx,Gx,~] = build_ops(Grid.x);
    [Dy,Gy,~] = build_ops(Grid.y);
    [Dp,Gp,~] = build_ops(Grid.p);

    %form divergence mats
    Dxx = Dx(:,1:Grid.x.Nfx); Dxy = Dx(:,(Grid.x.Nfx+1):Grid.x.Nf);
    Dyx = Dy(:,1:Grid.y.Nfx); Dyy = Dy(:,(Grid.y.Nfx+1):Grid.y.Nf);

    %form gradient mats
    Gxx = Gx(1:Grid.x.Nfx,:); Gxy = Gx((Grid.x.Nfx+1):Grid.x.Nf,:);
    Gyx = Gy(1:Grid.y.Nfx,:); Gyy = Gy((Grid.y.Nfx+1):Grid.y.Nf,:);
    
    %check if symmetric without natural boundary conditions
    %Gxx = Dxx'; Gxy = Dxy'; Gyx = Dyx'; Gyy = Dyy';
    
    %form eDot
    Edot = [Gxx,zeros(Grid.x.Nfx,Grid.y.N);
        zeros(Grid.y.Nfy,Grid.x.N),Gyy;
        .5*Gxy,.5*Gyx];


    D = [Dxx,zeros(Grid.x.N,Grid.y.Nfy),Dxy;
        zeros(Grid.y.N,Grid.x.Nfx),Dyy,Dyx];
    
    I = speye(Grid.N);




















