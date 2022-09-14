function u = get_value(x,geo,pd_u,pd_l,mu,open_or_closed )

   % this function gets the value of the field for the bump problem
    delta = 0.005;

    u = 0;

    Ns = 0;
    
    Nt = geo.get_N_total;
    
    xp = [];
    
    for p=1:geo.n_parts
        
        prt = geo.get_part(p);
        xp = [xp;prt.x];
        dw = [dw;prt.dw];
        
    end
    
    if strcmp (open_or_closed,'closed')
        isin = inpolygon(x(:,1),x(:,2),xp(:,1),xp(:,2));
    else strcmp (open_or_closed,'open')
        isin = inpolygon(x(:,1),x(:,2),[xp(:,1);xp(end,1);xp(1,1);xp(1,1)],[xp(:,2);-1e5;-1e5;xp(1,2)]);
    end
        
    %% only for the square
    
    
    for p=1:geo.n_parts
        
        prt = geo.get_part(p);
        
        [dist,xc,wc,rep] = prt.dist(x);
        
        if ~isin %rep==1
            
            k = pd_u.k;
            coef =1;
            
        else
            
            k = pd_l.k;
            coef =-1;
        end
        
        Np = prt.Np;
        
        if dist>delta

            D = DL(k,x,prt);
            
            S = SL(k,x,prt);
            
            u = u + coef*[D -S]*[prt.win.*mu(Ns+(1:Np));prt.win.*mu(Nt+Ns+(1:Np))];
            
        else

            ub = mu(Ns + 1:Np);

            ub = [ub(1);ub;ub(end)];

            w = [0;prt.w;2*pi];

            u = spline(w,ub,wc);   

            return

        end

        Ns = Ns  + prt.Np;

    end

end

%-------------------------------------------------------------------------%

function DL = DL(k,x,prt_i)
    
    Npe = size(x,1);
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
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);
    
    dW = repmat(prt_i.dw.',Npe,1);
    
    DL = (1i*k/4*besselh(1,k*r)./r.*((x_1-xp_1).*dxp_2-(x_2-xp_2).*dxp_1)*pi/prt_i.N).*dW;
    
    
end


function SL = SL(k,x,prt_i)
    
    Npi = prt_i.Np;
    Npe = size(x,1);

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

    tau = sqrt(dxp_1.^2+dxp_2.^2);
    
    dW = repmat(prt_i.dw.',Npe,1);
    
    SL = (1i/4.*besselh(0,k*r) .* tau * pi/prt_i.N).*dW;
    
    

end