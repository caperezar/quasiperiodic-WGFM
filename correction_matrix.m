function [M_c,Lp,Lm,dLp,dLm,proj_mod] = correction_matrix(pd,geo,H,N)

L = pd.L;
tol = 1e-10;

%% Matrix construction:
Nx  = round(N*L);
% Nx =100;
x0  = -L/2 +(0:Nx-1)*L/Nx;   
[dEp,Ep] = evaluate_matrix_horizontal(x0,H,pd,geo);
[dEm,Em] = evaluate_matrix_horizontal(x0,-H,pd,geo);

M_c = zeros(2*(geo.obstacles(1).Np*geo.n_o+geo.segments(1).Np));
r_1=[];
nr_1  = []; 
for n=1:geo.n_o
    r_1  = [r_1;geo.obstacles(n).x];
    nr_1 = [nr_1;geo.obstacles(n).normals];
end
    
np_2 = size(geo.segments(1).x,1);
n_reg_modes = numel(pd.reg_modes);

Ntot = 2*(geo.n_o*geo.obstacles(1).Np+geo.segments(1).Np);
Lp  = zeros(n_reg_modes,Ntot);
Lm  = zeros(n_reg_modes,Ntot);

flag = false;

for j=1:n_reg_modes 
    
    alpha_n = pd.reg_alpha_n(j);
    beta_n  = pd.reg_beta_n(j);

    u_p  = @(x) exp(1i*alpha_n*x(:,1)+1i*beta_n*x(:,2));
    du_p = @(x,nr) 1i*(alpha_n.*nr(:,1)+beta_n.*nr(:,2)).*u_p(x);

    u_m  = @(x) exp(1i*alpha_n*x(:,1)-1i*beta_n*x(:,2));
    du_m = @(x,nr) 1i*(alpha_n.*nr(:,1)-beta_n.*nr(:,2)).*u_m(x);

    Du_p  = @(x) 1i*exp(1i*alpha_n*x(:,1)+1i*beta_n*x(:,2)).*x(:,2);
    Ddu_p = @(x,nr) (-alpha_n * nr(:,1).*x(:,2) + nr(:,2).*(1i-beta_n.*x(:,2)))...
        .*exp(1i*alpha_n*x(:,1)+1i*beta_n*x(:,2));

    Du_m  = @(x) -1i*exp(1i*alpha_n*x(:,1)-1i*beta_n*x(:,2)).*x(:,2);
    Ddu_m = @(x,nr) (alpha_n * nr(:,1).*x(:,2) - nr(:,2).*(1i+beta_n.*x(:,2)))...
        .*exp(1i*alpha_n*x(:,1)-1i*beta_n*x(:,2));


    Phi_m = [u_p(r_1);du_p(r_1,nr_1);zeros(np_2,1);zeros(np_2,1)];
    Phi_p = [u_m(r_1);du_m(r_1,nr_1);zeros(np_2,1);zeros(np_2,1)];

    
    proj_mod = (1/Nx)*exp(-1i*alpha_n*x0);
    
    Lm(j,:) = exp(1i*beta_n*H)*proj_mod*(dEm+1i*beta_n*Em);
    Lp(j,:) = exp(1i*beta_n*H)*proj_mod*(dEp-1i*beta_n*Ep);

    if abs(beta_n)<=tol

        flag = true;

        dLm = (1/Nx)*exp(-1i*alpha_n*x0)*(1i*Em);
        dLp = (1/Nx)*exp(-1i*alpha_n*x0)*(-1i*Ep);
        
        DPhi_m = [Du_p(r_1);Ddu_p(r_1,nr_1);zeros(np_2,1);zeros(np_2,1)];
        DPhi_p = [Du_m(r_1);Ddu_m(r_1,nr_1);zeros(np_2,1);zeros(np_2,1)];

        M_c = M_c + exp(1i*beta_n*H)/(2i)*...
            (DPhi_m*Lm(j,:)-DPhi_p*proj_mod*Lp(j,:)+Phi_m*dLm -Phi_p*dLp);

    else
                  
        M_c = M_c + 1/(2i*beta_n)*(Phi_m*Lm(j,:)-Phi_p*Lp(j,:));
    end
       
end

if ~flag
    dLm = [];
    dLp = [];
end

end
function [dG_H,G_H] = evaluate_matrix_horizontal(x,H,pd,geo)
% Inputs: x, H, pd, and geo
%         the evaluation points have coordinates (x(i),H)

L = pd.L;
nu = pd.nu;
geo_large = geometry;
geo_large = geo_large.cctor(3*geo.n_o,geo.obstacles(1).name,2,geo.segments(1).name);

for n=1:geo.n_o

    geo_large.obstacles(n)           = geo.obstacles(n);
    geo_large.obstacles(n+geo.n_o)   = geo.obstacles(n);
    geo_large.obstacles(n+2*geo.n_o) = geo.obstacles(n);

    geo_large.obstacles(n).x(:,1)    = geo.obstacles(n).x(:,1)-L;
    geo_large.obstacles(n+2*geo.n_o).x(:,1) = geo.obstacles(n).x(:,1)+L;

end

geo_large.segments(1)  = geo.segments(1);
geo_large.segments(2)  = geo.segments(2);

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

A0 = zeros(Nx,2*geo.obstacles(1).Np*geo.n_o);
B0 = zeros(Nx,2*geo.obstacles(1).Np*geo.n_o);

for p=1:3
    D =  [];
    S =  [];

    dD =  [];
    dS =  [];
    for o=1:geo.n_o    

        o_ind = (p-1)*geo.n_o+o;

        pot_ext = pot_ext.cctor_local(pd.k(1),pts,geo_large.obstacles(o_ind));
        grad_pot_ext = grad_pot_ext.cctor_local(pd.k(1),pts,geo_large.obstacles(o_ind));
    
        dD  = [dD  gm^(p-2)*grad_pot_ext.DL(Nx+1:end,:)];
        dS =  [dS -nu*gm^(p-2)*grad_pot_ext.SL(Nx+1:end,:)];

        D  = [D  gm^(p-2)*pot_ext.DL];
        S  = [S -nu*gm^(p-2)*pot_ext.SL];

    end
    A0= A0 + [dD dS]; 
    B0= B0 + [D S];        
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
G_H  = [B0 B1];

end
