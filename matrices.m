classdef matrices
    properties
        
        SL % Single layer
        
        DL % Double layer
        
        DS % Normal derivative of single layer
        
        H % Hypersingular operators
        
        WD % Windowing operator (this is a diagonal matrix)  
        
        JCB
        
    end
    
    methods
        
%--------------------------------------------------------------------------------------------%    

        function obj = cctor_local(obj,k,geo,t_evl,t_int)
            
            prt_evl = geo.get_part(t_evl);
            
            prt_int  = geo.get_part(t_int);
            
            if (t_int==t_evl) 
                
                [obj.SL,obj.DL,obj.DS,obj.H] = self_interaction_matrix(k,prt_evl);
                
            else
                
                [obj.SL,obj.DL,obj.DS,obj.H] = far_interaction_matrix(k,prt_evl,prt_int);
                
            end
            
            obj.JCB = prt_int.dw;
            
            if geo.tags(t_int,1) ==2
                
                obj.WD = diag(prt_int.win);
                
            else
                
                obj.WD = eye(prt_int.Np);
                
            end
    
        end
        
%--------------------------------------------------------------------------------------------%

        function obj = cctor_self(obj,k,geo)
            
            
            for t_evl = 1:geo.n_parts
                
                prt_evl = geo.get_part(t_evl);
                
                for t_int = 1:geo.n_parts
                    
                    if (t_evl==t_int) 
                        
                        [S,D,DSL,HS] = self_interaction_matrix(k,prt_evl);
                        
                    else
                        
                        prt_int  = geo.get_part(t_int);
                        
                        [S,D,DSL,HS] = far_interaction_matrix(k,prt_evl,prt_int);
                        
                    end
                    
                    if (t_int == 1)

                        AS = S;
                        AD = D; 
                        ADS = DSL; 
                        AH = HS; 
                        
                    else

                        AS = [AS S];
                        AD = [AD D];
                        ADS = [ADS DSL];
                        AH = [AH HS];
                        
                    end

                end
                
                
                if (t_evl == 1)

                    BS = AS;
                    BD = AD;
                    BDS = ADS;
                    BH = AH;

                    BJ = prt_evl.dw;
                    
                    if geo.tags(t_evl,1) ==2 
                        
                        BW = prt_evl.win;
                        
                    else
                        
                        BW = ones(prt_evl.Np,1);
                        
                    end
                    
                
                else

                    BS = [BS;AS];
                    BD = [BD;AD];
                    BDS = [BDS;ADS];
                    BH = [BH;AH];
                    
                    BJ = [BJ;prt_evl.dw];
                                        
                    if geo.tags(t_evl,1) ==2 
                        
                        BW = [BW;prt_evl.win];
                        
                    else
                        
                        BW = [BW;ones(prt_evl.Np,1)];
                        
                    end
                                    
                end
                
            end
            
            obj.SL = BS;
            
            obj.DL = BD;
            
            obj.DS = BDS;
            
            obj.H = BH;
            
            obj.WD= diag(BW);
            
            obj.JCB = BJ;
           
        end
        

%-------------------------------------------------------------------------%
        function obj = cctor_far(obj,k,geo_evl,geo_int)
            
            for t_evl = 1:geo_evl.n_parts
                
                prt_evl = geo_evl.get_part(t_evl);
                
                for t_int = 1:geo_int.n_parts
                    
                    prt_int  = geo_int.get_part(t_int);
                        
                    [S,D,DSL,HS] = far_interaction_matrix(k,prt_evl,prt_int);
                    
                    if (t_int == 1)

                        AS = S;
                        AD = D; 
                        ADS = DSL; 
                        AH = HS; 
                        
                    else

                        AS = [AS S];
                        AD = [AD D];
                        ADS = [ADS DSL];
                        AH = [AH HS];
                        
                    end

                end
                
                
                if (t_evl == 1)

                    BS = AS;
                    BD = AD;
                    BDS = ADS;
                    BH = AH;
                    
                else

                    BS = [BS;AS];
                    BD = [BD;AD];
                    BDS = [BDS;ADS];
                    BH = [BH;AH];
                    
                end
                
            end
            
            obj.SL = BS;
            
            obj.DL = BD;
            
            obj.DS = BDS;
            
            obj.H = BH;
                       
        end
        

%-------------------------------------------------------------------------%
        function M = get_local_matrix(obj,geo,p_evl,p_int,type)
            
            Npi = geo.get_Np(p_int);
            Npe = geo.get_Np(p_evl);
            
            Npi_s = 0;
            
            if (p_int>1)
                
                for p = 1:p_int-1
                    
                    Npi_s = Npi_s + geo.get_Np(p);   
                    
                end
                
            end
            
            Npe_s = 0;
            
            if (p_evl>1)
                
                for p = 1:p_evl-1
                    
                    Npe_s= Npe_s + geo.get_Np(p);
                    
                end
                
            end
            
            
            if strcmp(type,'SL')
                
                M = obj.SL(Npe_s+(1:Npe),Npi_s+(1:Npi));
                
            elseif strcmp(type,'DL')
                
                M = obj.DL(Npe_s+(1:Npe),Npi_s+(1:Npi));
                
            elseif strcmp(type,'DS')
                
                M = obj.DS(Npe_s+(1:Npe),Npi_s+(1:Npi));
                
            elseif strcmp(type,'H')
                
                M = obj.H(Npe_s+(1:Npe),Npi_s+(1:Npi));
                
            end
            
        end
%-------------------------------------------------------------------------%    

        function obj = interaction_matrices(obj,k,prt_evl,prt_int)

            obj.SL = SL_far(k,prt_evl,prt_int);
            obj.DL = DL_far(k,prt_evl,prt_int);
            obj.DS = DS_far(k,prt_evl,prt_int);
            obj.H  = HS_far(k,prt_evl,prt_int);

        end

    end

end

function [S,D,DS,H] = self_interaction_matrix(k,prt)

    N = prt.N;

    R = nystrom_weights(N);
    
    S = SL_self(k,prt,R);

    D = DL_self(k,prt,R);
    
    DS = DS_self(k,prt,R);

    %if (strcmp(prt.name,'obstacle'))
    %    H = HS_self(k,prt,R);        
    %else
        H = NS_self(k,prt,R);
    %end

end

% ----------------------------------------------------------------------- %

function [S,D,DSL,H] = far_interaction_matrix(k,prt_evl,prt_int)

    S   = SL_far(k,prt_evl,prt_int);
    D   = DL_far(k,prt_evl,prt_int);
    DSL = DS_far(k,prt_evl,prt_int);
    H   = HS_far(k,prt_evl,prt_int);
    
end

%-------------------------------------------------------------------------------------%
function R = nystrom_weights(N)
    
    t = pi/N*(0:2*N-1);
    
    if N<=500
        tpmintj = repmat(t,2*N,1)-repmat(t.',1,2*N);

        temp_vec = zeros([1,1,N-1]);

        temp_vec(1,1,:)=(1:N-1);

        R=-2*pi/N*sum(cos(repmat(tpmintj,[1,1,N-1]).*repmat(temp_vec,[2*N,2*N,1]))./...
            repmat(temp_vec,[2*N,2*N,1]),3)-pi/(N^2)*cos(N*tpmintj);
    else

        R=zeros(2*N);
        for p=1:2*N
            temp=zeros(1,2*N);
            for j=1:2*N
                temp(j)=-2*pi/N*(sum((1./(1:N-1)).*cos((1:N-1)*(t(p)-t(j)))))-pi/(N^2)*cos(N*(t(p)-t(j)));
            end
            R(p,:)=temp;
        
        end

    end
end

%-------------------------------------------------------------------------%
%                       Single layer
%-------------------------------------------------------------------------%

function SL = SL_self(k,prt,R)
    
    t = prt.s;

    Np = prt.Np;

    x1 = prt.x(:,1);
    x2 = prt.x(:,2);

    dx1 = prt.dx(:,1);
    dx2 = prt.dx(:,2);

    xp_1 = repmat(x1.',Np,1);
    xp_2 = repmat(x2.',Np,1);

    x_1  = repmat(x1,1,Np);
    x_2   = repmat(x2,1,Np);

    dxp_1 = repmat(dx1.',Np,1);
    dxp_2 = repmat(dx2.',Np,1);
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);

    tau =sqrt(dxp_1.^2+dxp_2.^2);

    SL1 = (-1/(4*pi)).*besselj(0,k*r).*tau;

    SL2 = 1i/4.*besselh(0,k*r).*tau-SL1.*log(4*(sin((repmat(t,1,Np)-repmat(t.',Np,1))/2)).^2);

    dia = ~~eye(Np); 
    
    %size(tau)
    %size((1i/4+psi(1)/2/pi-1/4/pi*log(k^2/4*tau(dia).^2)).* tau(dia))   
    %size(2*log(prt.dw).*SL1(dia))
     
    SL2(dia)=(1i/4+psi(1)/2/pi-1/4/pi*log(k^2/4*tau(dia).^2)).* tau(dia)...
            +  2*log(prt.dw).*SL1(dia);
        
    dW = repmat(prt.dw.',Np,1);
    
    SL = (R(1:Np,1:Np).*SL1 + pi/prt.N*SL2).*dW;

end

%-------------------------------------------------------------------------%

function SL = SL_far(k,prt_e,prt_i)
    
    Npi = prt_i.Np;
    Npe = prt_e.Np;

    x1 = prt_e.x(:,1);
    x2 = prt_e.x(:,2);

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

%-------------------------------------------------------------------------%
%                             Double layer 
%-------------------------------------------------------------------------%

function DL = DL_self(k,prt,R)
    
    t = prt.s;

    Np = prt.Np;

    x1 = prt.x(:,1);
    x2 = prt.x(:,2);

    dx1 = prt.dx(:,1);
    dx2 = prt.dx(:,2);

    d2x1 = prt.d2x(:,1);
    d2x2 = prt.d2x(:,2);

    xp_1 = repmat(x1.',Np,1);
    xp_2 = repmat(x2.',Np,1);
    
    x_1  = repmat(x1,1,Np);
    x_2  = repmat(x2,1,Np);

    dxp_1 = repmat(dx1.',Np,1);
    dxp_2 = repmat(dx2.',Np,1);
    
    d2xp_1 = repmat(d2x1.',Np,1);
    d2xp_2 = repmat(d2x2.',Np,1);
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);

    tau =sqrt(dxp_1.^2+dxp_2.^2);

    dia = ~~eye(Np);
    
    DL1 = -k/(4*pi)*besselj(1,k*r)./r.*(dxp_2.*(x_1-xp_1)-dxp_1.*(x_2-xp_2));
    
    DL1(dia) = -k^2/(8*pi).*(dxp_2(dia).*(x_1(dia)-xp_1(dia))-dxp_1(dia).*(x_2(dia)-xp_2(dia)));
        
    
    DL2 = 1i*k/4*besselh(1,k*r)./r.*((x_1-xp_1).*dxp_2 - (x_2-xp_2).*dxp_1) ...
        -DL1.*log(4*(sin((repmat(t,1,Np)-repmat(t.',Np,1))/2)).^2);
        
    
        
    DL2(dia) = -1/(4*pi)*(dxp_1(dia).*d2xp_2(dia)-dxp_2(dia).*d2xp_1(dia))./tau(dia).^2 ...
            +2*log(prt.dw).*DL1(dia);
        
    dW = repmat(prt.dw.',Np,1);
    
    DL = (R(1:Np,1:Np).*DL1 + pi/prt.N*DL2).*dW;
    
end
%-------------------------------------------------------------------------%

function DL = DL_far(k,prt_e,prt_i)
    
    Npe = prt_e.Np;
    Npi = prt_i.Np;
    
    % Evaluation point
    
    x1 = prt_e.x(:,1);
    x2 = prt_e.x(:,2);

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

%-------------------------------------------------------------------------%
%                   Derivative of Single Layer
%-------------------------------------------------------------------------%

function DS = DS_self(k,prt,R)
    
    t = prt.s;

    Np = prt.Np;

    x1 = prt.x(:,1);
    x2 = prt.x(:,2);

    dx1 = prt.dx(:,1);
    dx2 = prt.dx(:,2);

    d2x1 = prt.d2x(:,1);
    d2x2 = prt.d2x(:,2);

    x_1  = repmat(x1,1,Np);
    x_2  = repmat(x2,1,Np);

    xp_1 = repmat(x1.',Np,1);
    xp_2 = repmat(x2.',Np,1);
    
    dx_1 = repmat(dx1,1,Np);
    dx_2 = repmat(dx2,1,Np);
    
    dxp_1 = repmat(dx1.',Np,1);
    dxp_2 = repmat(dx2.',Np,1);
    
    d2xp_1 = repmat(d2x1.',Np,1);
    d2xp_2 = repmat(d2x2.',Np,1);
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);

    taup = sqrt(dxp_1.^2 + dxp_2.^2);
    tau  = sqrt(dx_1.^2  + dx_2.^2);
    
    dia = ~~eye(Np);
    
    DS1 = k/(4*pi)*besselj(1,k*r)./r.*(dx_2.*(x_1-xp_1)-dx_1.*(x_2-xp_2)).*taup./tau;
    
    DS1(dia) = k^2/(8*pi).*(dx_2(dia).*(x_1(dia)-xp_1(dia))-dx_1(dia).*(x_2(dia)-xp_2(dia)));

    DS2 = 1i*k/4*besselh(1,k*r)./r.*((xp_1-x_1).*dx_2 - (xp_2-x_2).*dx_1).*taup./tau ...
        -DS1.*log(4*(sin((repmat(t,1,Np)-repmat(t.',Np,1))/2)).^2);


        
    DS2(dia) = -1/(4*pi)*(dxp_1(dia).*d2xp_2(dia)-dxp_2(dia).*d2xp_1(dia))./tau(dia).^2 ...
            +2*log(prt.dw).*DS1(dia);
    
    dW = repmat(prt.dw.',Np,1);
    
    DS = (R(1:Np,1:Np).*DS1 + pi/prt.N*DS2).*dW;
    
    
end
%-------------------------------------------------------------------------%
function DS = DS_far(k,prt_e,prt_i)

    Npe = prt_e.Np;
    Npi = prt_i.Np;
    
    % Evaluation point
    
    x1 = prt_e.x(:,1);
    x2 = prt_e.x(:,2);

    dx1 = prt_e.dx(:,1);
    dx2 = prt_e.dx(:,2);
    
    x_1  = repmat(x1,1,Npi);
    x_2  = repmat(x2,1,Npi);
    
    dx_1 = repmat(dx1,1,Npi);
    dx_2 = repmat(dx2,1,Npi);
    
    tau  = sqrt(dx_1.^2  + dx_2.^2);
    
    % Integration point
    
    xp1 = prt_i.x(:,1);
    xp2 = prt_i.x(:,2);
    
    dxp1 = prt_i.dx(:,1);
    dxp2 = prt_i.dx(:,2);
    xp_1 = repmat(xp1.',Npe,1);
    xp_2 = repmat(xp2.',Npe,1);
    
    dxp_1 = repmat(dxp1.',Npe,1);
    dxp_2 = repmat(dxp2.',Npe,1);
    
    taup  = sqrt(dxp_1.^2  + dxp_2.^2);
    
    %%%
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);
    
    dW = repmat(prt_i.dw.',Npe,1);
    
    DS = (1i*k/4*besselh(1,k*r)./r.*((xp_1-x_1).*dx_2 - (xp_2-x_2).*dx_1).*taup./tau * pi/prt_i.N).*dW;
    
end
%-------------------------------------------------------------------------%
%                           Hypersingular
%-------------------------------------------------------------------------%
function HS = HS_self(k,prt,R)    

    Np = prt.Np;

    I = eye(Np,Np);
    
    Der = FFTdiff(I,2*pi);
    
    S = SL_self(k,prt,R);
    
    HS = k^2*S.* (prt.normals * prt.normals');
    
    HS =  HS+diag(1./prt.tau)*Der*S*diag(1./prt.tau)*Der;
    
end

%-------------------------------------------------------------------------%

function HS = HS_far(k,prt_e,prt_i)

    Npe = prt_e.Np;
    Npi = prt_i.Np;
    
    % Evaluation point
    
    x1 = prt_e.x(:,1);
    x2 = prt_e.x(:,2);
    
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
    
    H1 = besselh(1,k*r);
    H0 = besselh(0,k*r);

    A11 = (1i/4)*k^2.*(x_1-xp_1).^2.*H0./r.^2 ...
        +(1i/4)*k*(r.^2-2*(x_1-xp_1).^2).*H1./r.^3;

    A12 = -(1i/2)*k*(x_1-xp_1).*(x_2-xp_2).*H1./r.^3 ...
        +(1i/4)*k^2*(x_1-xp_1).*(x_2-xp_2).*H0./r.^2;

    A21 = -(1i/2)*k*(x_1-xp_1).*(x_2-xp_2).*H1./r.^3 ...
        +(1i/4)*k^2*(x_1-xp_1).*(x_2-xp_2).*H0./r.^2;

    A22 = (1i/4)*k^2*(x_2-xp_2).^2.*H0./r.^2 ...
        +(1i/4)*k*(r.^2-2*(x_2-xp_2).^2).*H1./r.^3;
    
    n_x_1 = repmat(prt_e.normals(:,1),1,Npi);
    n_x_2 = repmat(prt_e.normals(:,2),1,Npi);
    
    n_xp_1 = dxp_2;
    n_xp_2 = -dxp_1;
    
    dW = repmat(prt_i.dw.',Npe,1);
    
    HS = (n_x_1.*(A11.*n_xp_1 + A12.*n_xp_2) + n_x_2.*(A21.*n_xp_1 + A22.*n_xp_2)).*dW*pi/prt_i.N;
    
    
    
    
%     r  = norm(x-xp);
%     H1 = besselh(1,k*r);
%     H0 = besselh(0,k*r);
% 
%     A(1,1) = (1i/4)*k^2*(x(1)-xp(1))^2*H0/r^2 ...
%     +(1i/4)*k*(r^2-2*(x(1)-xp(1))^2)*H1/r^3;
% 
%     A(1,2) = -(1i/2)*k*(x(1)-xp(1))*(x(2)-xp(2))*H1/r^3 ...
%     +(1i/4)*k^2*(x(1)-xp(1))*(x(2)-xp(2))*H0/r^2;
% 
%     A(2,1) = -(1i/2)*k*(x(1)-xp(1))*(x(2)-xp(2))*H1/r^3 ...
%     +(1i/4)*k^2*(x(1)-xp(1))*(x(2)-xp(2))*H0/r^2;
% 
%     A(2,2) = (1i/4)*k^2*(x(2)-xp(2))^2*H0/r^2 ...
%     +(1i/4)*k*(r^2-2*(x(2)-xp(2))^2)*H1/r^3;
%     
%     n_x = [dx(2) -dx(1)]/norm(dx);
%     n_xp = [dxp(2) -dxp(1)];
%     
%     HS  = n_x*A*n_xp';

    
    
    
end
%-------------------------------------------------------------------------%
function NS = NS_self(k,prt,R)
% This is the single layer operator without the singular kernel
% that is used in the solution of transmission problems

    t = prt.s;

    Np = prt.Np;

    x1 = prt.x(:,1);
    x2 = prt.x(:,2);

    dx1 = prt.dx(:,1);
    dx2 = prt.dx(:,2);

    d2x1 = prt.d2x(:,1);
    d2x2 = prt.d2x(:,2);

    d3x1 = prt.d3x(:,1);
    d3x2 = prt.d3x(:,2);
    
    % Evaluation point
    
    x_1  = repmat(x1,1,Np);
    x_2  = repmat(x2,1,Np);

    dx_1 = repmat(dx1,1,Np);
    dx_2 = repmat(dx2,1,Np);
    
    tau  = sqrt(dx_1.^2  + dx_2.^2);
    
    % Integration point
    
    xp_1 = repmat(x1.',Np,1);
    xp_2 = repmat(x2.',Np,1);
    
    
    dxp_1 = repmat(dx1.',Np,1);
    dxp_2 = repmat(dx2.',Np,1);
    
    %
    
    r = sqrt((x_1-xp_1).^2+(x_2-xp_2).^2);

    Nt = ((x_1-xp_1).*dxp_1+(x_2-xp_2).*dxp_2).*((x_1-xp_1).*dx_1+(x_2-xp_2).*dx_2)./r.^2;

    dia = ~~eye(Np);
    
    % NS1
    
    J0 = besselj(0,k*r);
    J1 = besselj(1,k*r);
    
    NS1 = -1/2/pi*Nt.*(k^2*J0-2*k*J1./r)-k/2/pi*(dx_1.*dxp_1+dx_2.*dxp_2)./r.*J1;
    
    NS1(dia) = -k^2/4/pi*prt.tau.^2;
    
    % NS
    
    H0 = besselh(0,k*r);
    H1 = besselh(1,k*r);
    
    NS = 1i/2*Nt .* ( k^2 * H0 - 2 * k * H1./r ) ...
            +1i*k/2*(dx_1.*dxp_1 + dx_2.*dxp_2)./r.*H1 ...
            +1/4/pi./(sin((repmat(t,1,Np)-repmat(t.',Np,1))/2)).^2;
    
        
    % NS2
    
    NS2 = NS-NS1.*log(4*(sin((repmat(t,1,Np)-repmat(t.',Np,1))/2)).^2);
    
    C = -psi(1);

    NS2(dia) = k^2/4/pi*(pi*1i-1-2*C-2*log(k*prt.tau/2)).*prt.tau.^2 +1/12/pi ...            
            +1/2/pi*(dx1.*d2x1+dx2.*d2x2).^2./prt.tau.^4 ...
            -1/4/pi*(d2x1.^2+d2x2.^2)./prt.tau.^2 ...
            -1/6/pi*(dx1.*d3x1+dx2.*d3x2)./prt.tau.^2;
        
    NS2(dia) = NS2(dia) +  2*log(prt.dw).*NS1(dia);
    
    %
    dW = repmat(prt.dw.',Np,1);
    
    NS = -(R(1:Np,1:Np).*NS1 + pi/prt.N*NS2).*dW./2./tau;
    
    S = SL_self(k,prt,R);
    
    NS = k^2*S.*( prt.normals * prt.normals.') + NS;
    
end

