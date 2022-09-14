function [U,X,Y] = evaluate_field(x,y,pd,geo,dens,Lp,Lm,dLp,dLm,opt)
L = pd.L;
nu = pd.nu;

Nx = numel(x(:));
Ny = numel(y(:));

[X,Y] = meshgrid(x,y);

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];

u = zeros(Nx*Ny,1);

gm = pd.gm;

%% Exterior domain
ind_1   = isin(pts,geo,1);     
pts_1   = pts(ind_1,:);
pot_ext = potentials;
F1 =reshape(dens.f1,geo.obstacles(1).Np,geo.n_o);
F2 =reshape(dens.f2,geo.obstacles(1).Np,geo.n_o);

for o=1:geo.n_o    
    pot_ext = pot_ext.cctor_local(pd.k(1),pts_1,geo.obstacles(o));
    u(ind_1) = u(ind_1) +  (pot_ext.DL*F1(:,o) - nu*pot_ext.SL*F2(:,o));   
end


pot_ext = pot_ext.cctor_local(pd.k(1),pts_1,geo.segments(1));
w = geo.segments(1).win;

u(ind_1) = u(ind_1) + (pot_ext.DL*(w.*dens.f3) - pot_ext.SL*(w.*dens.f4));

pot_ext = pot_ext.cctor_local(pd.k(1),pts_1,geo.segments(2));
u(ind_1) = u(ind_1) +  (-pot_ext.DL*(w.*dens.f5) + pot_ext.SL*(w.*dens.f6));


% Correction:
f = [dens.f1;dens.f2;dens.f3;dens.f4];
lp = Lp*f;
lm = Lm*f;

if ~isempty(dLp)
    dlp = dLp*f;
    dlm = dLm*f;
end

for j=1:numel(pd.reg_beta_n)

    alpha_n = pd.reg_alpha_n(j);
    beta_n  = pd.reg_beta_n(j);
        
    u_p  = @(x) exp(1i*alpha_n*x(:,1)+1i*beta_n*x(:,2));
    u_m  = @(x) exp(1i*alpha_n*x(:,1)-1i*beta_n*x(:,2));
    Du_p  = @(x) 1i*exp(1i*alpha_n*x(:,1)+1i*beta_n*x(:,2)).*x(:,2);
    Du_m  = @(x) -1i*exp(1i*alpha_n*x(:,1)-1i*beta_n*x(:,2)).*x(:,2);
    
    if ~isempty(dLp)
        u(ind_1) = u(ind_1)+ 1/2i*(-dlm(j)*u_p(pts_1) + dlp(j)*u_m(pts_1)...
                                   -lm(j)*Du_p(pts_1) + lp(j)*Du_m(pts_1));
    else        
        u(ind_1) = u(ind_1)+ 1/(2i*beta_n)*(-lm(j)*u_p(pts_1)+lp(j)*u_m(pts_1));    
    end
end

if strcmp(opt,'t')
    
   u(ind_1) = u(ind_1) + exp(1i*pd.alpha*pts_1(:,1)+1i*pd.beta*pts_1(:,2));   
   
end

%% Interior domain

ind_2 = isin(pts,geo,2);     

pts_2 = pts(ind_2,:);

pot_int = potentials;

for o=1:geo.n_o    
    
    pot_int = pot_int.cctor_local(pd.k(2),pts_2,geo.obstacles(o));

    u(ind_2) =u(ind_2) + (-pot_int.DL*F1(:,o) + pot_int.SL*F2(:,o));

end


U = reshape(u,Ny,Nx);
 

end

% ----------------------------------------------------------------------- %

function in = isin(x,geo,dom)
    xv = geo.obstacles(1).x;
    in = inpolygon(x(:,1),x(:,2),xv(:,1),xv(:,2));
    for n=2:geo.n_o
        xv = geo.obstacles(n).x;
        in = or(in,inpolygon(x(:,1),x(:,2),xv(:,1),xv(:,2)));
    end
    
    if dom==1
        in = ~in;
    end

end

