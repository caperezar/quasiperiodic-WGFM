classdef potentials
    properties
        
        SL % Single layer
        DL % Double layer

    end
    
    methods

         function obj = cctor_local(obj,k,x,prt)
            
            obj.SL = SL_potential(k,x,prt);
            obj.DL = DL_potential(k,x,prt);
             
        end
        
        function obj = cctor(obj,k,x,geo)
            
            for t_int = 1:geo.n_parts
                    
                prt_int  = geo.get_part(t_int);

                S = SL_potential(k,x,prt_int);
                D = DL_potential(k,x,prt_int);
               
                if (t_int == 1)

                    AS = S;
                    AD = D; 

                else

                    AS = [AS S];
                    AD = [AD D];
                    
                end

            end

            obj.SL = AS;
            obj.DL = AD;
             
    end

%-------------------------------------------------------------------------%
    end

end


%%
function DL = DL_potential(k,x,prt_i)
    
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
    %tol  = 1e-4;
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);
    
    %r = r.* (r>tol) + tol*(r<=tol);
    
    dW = repmat(prt_i.dw.',Npe,1);
    
    DL = (1i*k/4*besselh(1,k*r)./r.*((x_1-xp_1).*dxp_2-(x_2-xp_2).*dxp_1)*pi/prt_i.N).*dW;

    %% Complex interpolation
%     tol  = 1e-4;
%     
%     ind = (r<=tol);
%     
%     h = 2*tol;
    
%     for n=1:size(ind,1)
%         
%         for m=1:size(ind,2)
%             
%             if ind(n,m)
%                 
%                 x_eval = x(n,:);
%                 x_int = prt.x(m,:);
%                 
%                 x1 = x_eval+1i*h/2*[1 1];
%                 r1 = sqrt((x1(1)-x_int(1)).^2+(x1(2)-x_int(2)).^2);
%                 A1 = besselh(1,k*r1)/r1;
%                 
%                 x2 = x_eval+1i*h/2*[-1 1];
%                 r2 = sqrt((x2(1)-x_int(1)).^2+(x2(2)-x_int(2)).^2);
%                 A2 = besselh(1,k*r2)/r2;
%                 
%                 x3 = x_eval+1i*h/2*[-1 -1];
%                 r3 = sqrt((x4(1)-x_int(1)).^2+(x3(2)-x_int(2)).^2);
%                 A3 = besselh(1,k*r3)/r3;
%                 
%                 x4 = x_eval+1i*h/2*[1 -1];
%                 r4 = sqrt((x4(1)-x_int(1)).^2+(x4(2)-x_int(2)).^2);
%                 A4 = besselh(1,k*r4)/r4;
%                 
%                 DL(n,m) = (1i/16.*(A1+A2+A3+A4)*((x_1(n,m)-xp_1(n,m))*dxp_2(n,m)-(x_2(n,m)-xp_2(n,m))*dxp_1(n,m))*pi/prt_i.N)*dW(n,m);
%                 
%             end                 
%             
%         end
        
    %end
    
end
%%
function SL = SL_potential(k,x,prt_i)
    
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
    
    tol  = 1e-4;
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);
    
    %r = r.* (r>tol) + tol*(r<=tol);

    tau = sqrt(dxp_1.^2+dxp_2.^2);
    
    dW = repmat(prt_i.dw.',Npe,1);
    
    SL = (1i/4.*besselh(0,k*r) .* tau * pi/prt_i.N).*dW;
    
    %% Complex interpolation
%     tol  = 1e-4;
%     
%     ind = (r<=tol);
%     
%     h = 2*tol;
%     
%     for n=1:size(ind,1)
%         
%         for m=1:size(ind,2)
%             
%             if ind(n,m)
%                 
%                 x_eval = x(n,:);
%                 x_int = prt.x(m,:);
%                 
%                 x1 = x_eval+1i*h/2*[1 1];
%                 r1 = sqrt((x1(1)-x_int(1)).^2+(x1(2)-x_int(2)).^2);
%                 A1 = besselh(0,k*r1);
%                 
%                 x2 = x_eval+1i*h/2*[-1 1];
%                 r2 = sqrt((x2(1)-x_int(1)).^2+(x2(2)-x_int(2)).^2);
%                 A2 = besselh(0,k*r2);
%                 
%                 x3 = x_eval+1i*h/2*[-1 -1];
%                 r3 = sqrt((x4(1)-x_int(1)).^2+(x3(2)-x_int(2)).^2);
%                 A3 = besselh(0,k*r3);
%                 
%                 x4 = x_eval+1i*h/2*[1 -1];
%                 r4 = sqrt((x4(1)-x_int(1)).^2+(x4(2)-x_int(2)).^2);
%                 A4 = besselh(0,k*r4);
%                 
%                 SL(n,m) = (1i/16.*(A1+A2+A3+A4)* tau(n,m) * pi/prt_i.N)*dW(n,m);
%                 
%             end                 
%             
%         end
%         
%     end
    
end
