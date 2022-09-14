classdef geometry
    properties
        obstacles
        segments
        n_o  % number of obstacles
        n_s  % number of segments
        tags % 1=obstacle,2=segment
        n_parts        

        top_seg
        bottom_seg
    end
    
    methods
%-------------------------------------------------------------------------%        
        function obj = cctor(obj,n_o,o_type,n_s,s_type)

            if nargin==3
                s_type = 'line_segment';
            end

            obj.n_o = n_o;
            obj.n_s = n_s;
                        
            obj.n_parts = n_o+n_s;
            obj.tags= ones(obj.n_parts,2);
            
            if n_o>0
                b(1,n_o) = eval(o_type);
                obj.obstacles =b;
                obj.tags(1:n_o,2) = (1:n_o)'; 
            end
            
            if n_s>0
                s(1,n_s) = eval(s_type);
                obj.segments =s;
                obj.tags(n_o+1:n_o+n_s,1) = 2*obj.tags(n_o+1:n_o+n_s,1);
                obj.tags(n_o+1:n_o+n_s,2) = (1:n_s)';
            end
            
        end
%-------------------------------------------------------------------------%        
        function  N = get_N_total(obj)
            N = 0;
            if obj.n_o>0
                for i=1:obj.n_o
                    N=N + 2*obj.obstacles(i).N;
                end
            end
            if obj.n_s>0
                for i=1:obj.n_s
                    N=N + 2*obj.segments(i).N-1;
                end
            end
            
        end
%-------------------------------------------------------------------------%        
        function prt = get_part(obj,t)
            % this function gets the t part of the geometry
            l = obj.tags(t,2);    
            if obj.tags(t,1) == 1 % it's an obstacle        
                prt = obj.obstacles(l);
            elseif obj.tags(t,1) ==2 % it's a line segment
                prt = obj.segments(l);
            end
        end
%-------------------------------------------------------------------------%
        function Np = get_Np(obj,t)
            % this function computes the number of points on the t part of
            % the geometry
            
            l = obj.tags(t,2);    
            if obj.tags(t,1) == 1 % it's an obstacleﬂﬂ
                Np = obj.obstacles(l).Np;
            elseif obj.tags(t,1) ==2 % it's a line segment
                Np = obj.segments(l).Np;
            end
        end
        
%-------------------------------------------------------------------------%
        function t_indices = get_indices(obj,t)
            % this function computes the indices of the points of the t
            % part of the goemetry
           
            n0 = 0;
            
            for p=1:t-1
                n0  = n0 + get_Np(obj,p);
            end
            
            ni = n0 + 1;
            
            ne = n0 + get_Np(obj,t);
            
            t_indices = (ni:ne)';
            
        end
%-------------------------------------------------------------------------%
        function win = get_window(obj)
            win = [];
            for t = 1:obj.n_parts
                prt = obj.get_part(t);
                win = [win;prt.win];                
            end
        end
%-------------------------------------------------------------------------%
function  plot(obj,nfig)
            figure(nfig); hold on
            for t = 1:obj.n_parts
                prt = obj.get_part(t);
                if strcmp(prt.name,'obstacle')
                    plot([prt.x(:,1);prt.x(1,1)],[prt.x(:,2);prt.x(1,2)]);
                else
                    plot(prt.x(:,1),prt.x(:,2));
                end
            end
            axis equal 
            axis tight
        end
%-------------------------------------------------------------------------%
    end
end
