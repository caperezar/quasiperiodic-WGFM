 classdef rhs 
    properties
        
        IF    % incident field
        dIF   % normal derivative of the incident field    
        
    end
    
    methods
%-------------------------------------------------------------------------%        
        function obj = cctor_local(obj,problem_data,geo,tag)
            prt = geo.get_part(tag);  
            [obj.IF,obj.dIF] = boundary_data(problem_data,prt); 
        end
%-------------------------------------------------------------------------%

        function obj = cctor_global(obj,problem_data,geo)
            
            for p = 1:geo.n_parts
                
                prt = geo.get_part(p);  
                
                [fi,dfi] = boundary_data(problem_data,prt);  
                
                if (p==1)
                    
                    FI  = fi;
                    dFI = dfi;

                else
                    
                    FI  = [FI;fi];
                    dFI = [dFI;dfi];
                    
                end

            end

            obj.IF  = FI;
            obj.dIF = dFI;

        end

    end
    
end

% ----------------------------------------------------------------------- %

function [fI,dfI] = boundary_data(problem_data,prt)

    Np = prt.Np;
    
    fI = zeros(Np,1);
    dfI = zeros(Np,1); 
    
    for i=1:Np
    
        x  = prt.x(i,:);
    
        n = prt.normals(i,:);

        [fI(i),du] = incident_field(x,problem_data); 
        
        dfI(i) = du(1)*n(1)+du(2)*n(2);
        
    end

    
end

