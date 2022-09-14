classdef line_segment
    properties
        x_start % x coordenate of the starting point
        x_end % x coordenate of the end point
        N % level of discretization (the total number of point on the bump is 
        p % degree of polinomial change of variable for Colton and Kress method
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
        function [x, dx,d2x,d3x] = param(obj,t)
            x_s = obj.x_start;
            x_e = obj.x_end;
        % Parametrization of the curve between 0:2*pi
            x = x_s+(x_e-x_s)*t/(2*pi);
            dx = (x_e-x_s)/(2*pi);
            d2x = [0 0];
            d3x = [0 0];
        end
        
%--------------------------------------------------------------------------------------------%        
        function obj = cctor(obj,N,x_s,x_e,percent)

%             obj.x_start = x_s;
% 
%             obj.x_end   = x_e;
% 
%             obj.N  = N;
%             
%                            
%             obj.Np  = 2*N-1;
%             obj.p  = 4;
%             obj.s = pi/N*(1:1:2*N-1)';
%             [obj.w,obj.dw] = W(obj.s,obj.p);
%             Nt = 2*N-1;
%                 
%                         
%             obj.x   = zeros(Nt,2);
%             obj.dx  = zeros(Nt,2);
%             obj.d2x = zeros(Nt,2);
%             obj.d3x = zeros(Nt,2);
%             
%             obj.win = zeros(Nt,1);  
%             obj.normals  = zeros(Nt,2); 
%             obj.tau =  zeros(Nt,1); 
%             
%             for i=1:Nt
% 
%                 [obj.x(i,:),obj.dx(i,:),obj.d2x(i,:),obj.d3x(i,:)] = obj.param(obj.w(i));
%                 
%                 x0 = 0.5*(x_s+x_e);
%                 A = norm(x_e-x0);
% 
%                 obj.win(i) =POU(norm(obj.x(i,:)-x0),A*percent,A);                                                                                                               
%                 
%                 obj.tau(i) = norm(obj.dx(i,:));
%                 
%                 obj.normals(i,:) = [obj.dx(i,2) -obj.dx(i,1)]/obj.tau(i);
%                 
%             end
%             
%             obj.name = 'line_segment'; 
                        
            
            obj.x_start = x_s;

            obj.x_end   = x_e;

            obj.N  = N;
                                       
            obj.Np  = 2*N;
            obj.s = pi/N*(0:1:2*N-1)';
            obj.w = obj.s;
            obj.dw = ones(size(obj.s));
            Nt = 2*N;
            
            
            obj.x   = zeros(Nt,2);
            obj.dx  = zeros(Nt,2);
            obj.d2x = zeros(Nt,2);
            obj.d3x = zeros(Nt,2);
            
            obj.win = zeros(Nt,1);  
            obj.normals  = zeros(Nt,2); 
            obj.tau =  zeros(Nt,1); 
            
            for i=1:Nt

                [obj.x(i,:),obj.dx(i,:),obj.d2x(i,:),obj.d3x(i,:)] = obj.param(obj.w(i));
                
                
                x0 = 0.5*(x_s+x_e);
                A = norm(x_e-x0);
                obj.win(i) =POU(norm(obj.x(i,:)-x0),A*percent,A);
                                                                                                      
                obj.tau(i) = norm(obj.dx(i,:));
                
                obj.normals(i,:) = [obj.dx(i,2) -obj.dx(i,1)]/obj.tau(i);
                
            end
            
            obj.name = 'line_segment'; 
            
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