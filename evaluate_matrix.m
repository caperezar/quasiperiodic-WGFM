function [dG_H,G_H] = evaluate_matrix(x,H,pd,geo)
% Inputs: x, H, pd, and geo
%         the evaluation points have coordinates (x(i),H)

L = pd.L;
geo_large = geometry;
geo_large = geo_large.cctor(3,2);

geo_large.obstacles(1) = geo.obstacles(1);
geo_large.obstacles(2) = geo.obstacles(1);
geo_large.obstacles(3) = geo.obstacles(1);
geo_large.segments(1)  = geo.segments(1);
geo_large.segments(2)  = geo.segments(2);

geo_large.obstacles(1).x(:,1) = geo.obstacles(1).x(:,1)-L;
geo_large.obstacles(3).x(:,1) = geo.obstacles(1).x(:,1)+L;

geo_large.segments(1).x(:,1)  = geo.segments(1).x(:,1)-L;
geo_large.segments(1).x_start(:,1)  = geo.segments(:,1).x_start(1)-L;
geo_large.segments(1).x_end(:,1)  = geo.segments(:,1).x_end(1)-L;


geo_large.segments(2).x(:,1)  = geo.segments(2).x(:,1)+L;
geo_large.segments(2).x_start(:,1)  = geo.segments(2).x_start(:,1)+L;
geo_large.segments(2).x_end(:,1)  = geo.segments(2).x_end(:,1)+L;


Nx = numel(x(:));

pts = [x' repmat(H,Nx,1)];

pot_ext = potentials;
grad_pot_ext = grad_potentials;

gm = pd.gm;

A0 = zeros(Nx,2*geo.obstacles.Np);
B0 = zeros(Nx,2*geo.obstacles.Np);
for o=1:3
    pot_ext = pot_ext.cctor_local(pd.k(1),pts,geo_large.obstacles(o));
    grad_pot_ext = grad_pot_ext.cctor_local(pd.k(1),pts,geo_large.obstacles(o));
    
    A0  = A0+gm^(o-2)*[grad_pot_ext.DL(Nx+1:end,:) -grad_pot_ext.SL(Nx+1:end,:)];
    B0  = B0+gm^(o-2)*[pot_ext.DL -pot_ext.SL];
    
end


pot_ext = pot_ext.cctor_local(pd.k(1),pts,geo_large.segments(1));
grad_pot_ext = grad_pot_ext.cctor_local(pd.k(1),pts,geo_large.segments(1));
W = diag(geo_large.segments(1).win);

A1  = gm^(-1)*[grad_pot_ext.DL(Nx+1:end,:)*W -grad_pot_ext.SL(Nx+1:end,:)*W];
B1 =  gm^(-1)*[pot_ext.DL*W -pot_ext.SL*W];

pot_ext = pot_ext.cctor_local(pd.k(1),pts,geo_large.segments(2));
grad_pot_ext = grad_pot_ext.cctor_local(pd.k(1),pts,geo_large.segments(2));

A1 = A1 +  gm^2*[-grad_pot_ext.DL(Nx+1:end,:)*W grad_pot_ext.SL(Nx+1:end,:)*W];
B1 = B1 +  gm^2*[-pot_ext.DL*W pot_ext.SL*W];

dG_H = [A0 A1];
G_H = [B0 B1];

end