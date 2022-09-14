function [Error,R,T] = energy_test(H,geo,pd,dens,Lp,Lm,dLp,dLm)

L = pd.L;

Nx = 100;h= L/Nx;

x0 = -L/2 +(0:Nx-1)*h;   

u = evaluate_field_super_cell(x0,[-H H],pd,geo,dens,Lp,Lm,dLp,dLm,'s');

beta = abs(pd.beta);
S = 0.0;
R = 0;
% T = 1;

for j=1:numel(pd.beta_n)

    alpha_n = pd.alpha_n(j);
    beta_n  = pd.beta_n(j);

    apn = h/L*exp(-1i*alpha_n*x0)*u(2,:).';
    Bp = exp(-1i*beta_n*H)*apn;
                    
    amn = h/L*exp(-1i*alpha_n*x0)*u(1,:).';
    Bm = exp(-1i*beta_n*H)*amn;
    
    if pd.prop_modes(j)==0
        B0 = Bm; 
    end
    S = S + beta_n/beta*(abs(Bp)^2+abs(Bm)^2);
    R = R + beta_n/beta*abs(Bp)^2;
%     T = T + beta_n/beta*abs(Bm)^2;

%     R = R + abs(Bp)^2;
%     T = T + abs(Bm)^2;
end

Error = abs(2*real(B0)+S);
% T = T+2*real(B0);
T=1+2*real(B0)+S-R;
% R = R/2/real(B0);
% T = T/2/real(B0);

