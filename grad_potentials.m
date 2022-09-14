classdef grad_potentials
    
    properties
        
        SL % Single layer
        
        DL % Double layer

    end
    
    methods

        function obj = cctor_local(obj,k,x,prt)
            
            obj.SL = DSL_potential(k,x,prt);
            obj.DL = DDL_potential(k,x,prt);
             
        end
        
        function obj = cctor(obj,k,x,geo)
            
            for t_int = 1:geo.n_parts
                    
                prt_int  = geo.get_part(t_int);

                DS = DSL_potential(k,x,prt_int);
                DD = DDL_potential(k,x,prt_int);
                 
                     
                if (t_int == 1)

                    AS = DS;
                    AD = DD; 

                else

                    AS = [AS DS];
                    AD = [AD DD];
                    
                end

            end

        obj.SL = AS;

        obj.DL = AD;
                    
    end

%-------------------------------------------------------------------------%
    end

end


%%
function DSL = DSL_potential(k,x,prt_i)
    
    Npi = prt_i.Np;
    Npe = numel(x(:,1));

    x1 = x(:,1);
    x2 = x(:,2);

    xp1 = prt_i.x(:,1);
    xp2 = prt_i.x(:,2);
    
    dxp1 = prt_i.dx(:,1);
    dxp2 = prt_i.dx(:,2);

    xp_1 = repmat(xp1.',Npe,1);
    xp_2 = repmat(xp2.',Npe,1);
    
    x_1  = repmat(x1,1,Npi);
    x_2  = repmat(x2,1,Npi);

    dxp_1 = repmat(dxp1.',Npe,1);
    dxp_2 = repmat(dxp2.',Npe,1);
    
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);
    
    %tol  = 1e-2;
    
    %r = r.* (r>tol) + tol*(r<=tol);

    tau = sqrt(dxp_1.^2+dxp_2.^2);
    
    dW = repmat(prt_i.dw.',Npe,1);
    
    A = 1i*k/4*besselh(1,k*r)./r;
    
    DSL = -A.*(x_1-xp_1).*tau.*dW*pi/prt_i.N;
    
    DSL = [DSL;-A.*(x_2-xp_2).*tau.*dW*pi/prt_i.N];
    
end

%%
function DDL = DDL_potential(k,x,prt_i)
    
   Npe = numel(x(:,1));
    Npi = prt_i.Np;
    
    % Evaluation point
    
    x1 = x(:,1);
    x2 = x(:,2);

    x_1  = repmat(x1,1,Npi);
    x_2  = repmat(x2,1,Npi);
    
    % Integration point
    
    xp1 = prt_i.x(:,1);
    xp2 = prt_i.x(:,2);
    
    dxp1 = prt_i.dx(:,1);
    dxp2 = prt_i.dx(:,2);
    xp_1 = repmat(xp1.',Npe,1);
    xp_2 = repmat(xp2.',Npe,1);
    
    dxp_1 = repmat(dxp1.',Npe,1);
    dxp_2 = repmat(dxp2.',Npe,1);
    
    %%%
    %tol  = 1e-2;
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);
    
    %r = r.* (r>tol) + tol*(r<=tol);
    
    dW = repmat(prt_i.dw.',Npe,1);%.* (r>tol);
    
    A = 1i*k/4*besselh(1,k*r)./r;
    
    dAdr = (1i/4)*k*(besselh(0, k*r).*k.*r-2*besselh(1, k*r))./r.^2;
    
    drdx1 = (x_1-xp_1)./r;
    drdx2 = (x_2-xp_2)./r;
    dAdx1 = dAdr.*drdx1;
    dAdx2 = dAdr.*drdx2;
    
    %DL = (1i*k/4*besselh(1,k*r)./r.*((x_1-xp_1).*dxp_2-(x_2-xp_2).*dxp_1)*pi/prt_i.N).*dW;
    
    %DDL = (A.*((x_1-xp_1).*dxp_2-(x_2-xp_2).*dxp_1)*pi/prt_i.N).*dW;
    DDL = (dAdx1.*((x_1-xp_1).*dxp_2-(x_2-xp_2).*dxp_1) + A.*dxp_2).*dW*pi/prt_i.N;
    DDL = [DDL; (dAdx2.*((x_1-xp_1).*dxp_2-(x_2-xp_2).*dxp_1) - A.*dxp_1).*dW*pi/prt_i.N];
    
    
end
    
    
%end