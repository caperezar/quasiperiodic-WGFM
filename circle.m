classdef circle
    properties
        loc % center of the obstacle
        rad % radiuous of the obstacle
        
        N % level of discretization (the total number of point on the bump is 2*N)
        s
        w
        dw
        
        x  
        dx 
        d2x
        d3x
        
        Np 
        win
        normals
        tau
        name
        
    end
   
    methods
  %--------------------------------------------------------------------------------------------%      
        function obj = cctor(obj,N,loc,rad)
            
            obj.N  = N;
            
            obj.Np  = 2*N;
            
            obj.loc = loc;

            obj.rad = rad;
    
            obj.s = pi/N*(0:1:2*N-1)';
            
            obj.w = obj.s;
            
            obj.dw = ones(size(obj.s));
            
            Nt = 2*N;
    
            obj.x  = zeros(Nt,2);
            obj.dx = zeros(Nt,2);
            obj.d2x = zeros(Nt,2);
            obj.d3x = zeros(Nt,2);
            
            obj.normals  = zeros(Nt,2); 
            obj.tau =  zeros(Nt,1); 
            obj.win =  ones(Nt,1); 
            
            for i=1:Nt
        
                [obj.x(i,:),obj.dx(i,:),obj.d2x(i,:),obj.d3x(i,:)] = obj.param(obj.s(i));
                obj.tau(i) = norm(obj.dx(i,:));
                
                obj.normals(i,:) = [obj.dx(i,2) -obj.dx(i,1)]/obj.tau(i);
                
            end
            
            obj.name = 'circle';
            
        end
 %--------------------------------------------------------------------------------------------%      
        function [x, dx,d2x,d3x] = param(obj,t)
        % Original parametrization of the curve between 0:2*pi
            %circle
            x  = obj.loc+obj.rad*[cos(t) sin(t)];
            dx = obj.rad*[-sin(t) cos(t)];
            d2x = -obj.rad*[cos(t) sin(t)];
            d3x = -obj.rad*[-sin(t) cos(t)];       

        end
        
%-------------------------------------------------------------------------%
        function [dist,xc,tc,rel_pos] = dist(obj,x)
            % Position of the point x with respect to the curve is
            % defined as follows:
            %
            %               x
            % Obove     <------< : rel_pos = 1
            %
            %
            % Below     <------< : rel_pos = -1
            %               x
            
            options = optimset ('TolX',1.0e-5); 
            
            f = @(t) norm(obj.param(t)-x);
            
            [tc,dist] = fminbnd(f,0,2*pi,options);
            
            [xc,dxc] = obj.param(tc);
            
            dxc = [dxc 0];
            d = [x-xc 0];
            
            A = cross(d,dxc);
            
            rel_pos = sign(A(3));
            
        end
        
%-------------------------------------------------------------------------%        
    end
    
end