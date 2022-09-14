function [M,W,E] = create_matrices_periodic(pd,geo,opt)
if nargin ==2
    opt = 'general';
end

gm = pd.gm;
nu = pd.nu;
k  = pd.k;

%% Obstacle-obstacle interactions
m1 = matrices;
m2 = matrices;

A11 = [];
A12 = [];
A21 = [];
A22 = [];

for n=1:geo.n_o
    A11_aux = [];
    A12_aux = [];
    A21_aux = [];
    A22_aux = [];
    for m=1:geo.n_o
        if m==n
            m1 = m1.cctor_local(k(1),geo,n,m);
            m2 = m2.cctor_local(k(2),geo,n,m);
            A11_aux = [A11_aux  m2.DL-m1.DL];
            A12_aux = [A12_aux -m2.SL+nu*m1.SL]; 
            A21_aux = [A21_aux  m2.H-m1.H];
            A22_aux = [A22_aux -m2.DS+nu*m1.DS];
        else
            m1 = m1.cctor_local(k(1),geo,n,m);
            A11_aux = [A11_aux  -m1.DL];
            A12_aux = [A12_aux nu*m1.SL]; 
            A21_aux = [A21_aux  -m1.H];
            A22_aux = [A22_aux +nu*m1.DS];
        end
    end
    A11 = [A11;A11_aux];
    A12 = [A12;A12_aux];
    A21 = [A21;A21_aux];
    A22 = [A22;A22_aux];
end
T11 = [A11 A12;A21 A22];

%% Obstacle-open curve interactions
m12 = matrices;
m13 = matrices;
A11 = [];
A12 = [];
A21 = [];
A22 = [];
for n=1:geo.n_o
    m12 = m12.cctor_local(k(1),geo,n,geo.n_o+1);
    m13 = m13.cctor_local(k(1),geo,n,geo.n_o+2);
    A11 = [A11;-m12.DL+gm*m13.DL];
    A12 = [A12; m12.SL-gm*m13.SL];  
    A21 = [A21;-m12.H+gm*m13.H];
    A22 = [A22; m12.DS-gm*m13.DS];   
end
T12 = [A11 A12;A21 A22];

%% Open curve-obstacle interactions
m21 = matrices;
m31 = matrices;
A11 = [];
A12 = [];
A21 = [];
A22 = [];
for n=1:geo.n_o
    m21 = m21.cctor_local(k(1),geo,geo.n_o+1,n);
    m31 = m31.cctor_local(k(1),geo,geo.n_o+2,n);
    A11 = [A11 -(m31.DL+gm*m21.DL)];
    A12 = [A12 nu*(m31.SL+gm*m21.SL)];        
    A21 = [A21 -(m31.H+gm*m21.H)];
    A22 = [A22 nu*(m31.DS+gm*m21.DS)];
    
end
T21 = [A11 A12;A21 A22];

%% Open curve-open curve interactions
m32 = matrices;
m32 = m32.cctor_local(k(1),geo,geo.n_o+2,geo.n_o+1);
if strcmp(opt,'general')
    m23 = matrices; 
    m23 = m23.cctor_local(k(1),geo,geo.n_o+1,geo.n_o+2);
    T22 = [-(m32.DL-gm^2*m23.DL)  (m32.SL-gm^2*m23.SL);
            -(m32.H-gm^2*m23.H)   (m32.DS-gm^2*m23.DS)];
elseif strcmp(opt,'line')
    T22 = [-(1+gm^2)*m32.DL  (1-gm^2)*m32.SL;
           -(1-gm^2)*m32.H   (1+gm^2)*m32.DS];
end
    
M = [T11 T12;T21 T22];

%% Window matrix
np_obs = geo.n_o*geo.obstacles(1).Np;
W11 = eye(2*np_obs);
W22 = diag([geo.segments(1).win;geo.segments(1).win]);
W   = blkdiag(W11,W22);

%% Invertible term
E11 = diag([ones(np_obs,1);0.5*(1+nu)*ones(np_obs,1)]);
E22 = gm * eye(2*geo.segments(1).Np);

E = blkdiag(E11,E22);








 
 
 