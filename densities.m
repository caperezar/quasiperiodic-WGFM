classdef densities    
    properties
        
     f1 % on the obstacle's boundary
     f2

     f3 % on the left wall of the unit cell 
     f4

     f5 % on the right wall of the unit cell
     f6 

    end
    
    methods
        
        function obj = cctor(obj,f,pd,geo)

            obs_indices = [];

            for n=1:geo.n_o            

                obs_indices = [obs_indices;geo.get_indices(n)];

            end

            seg_indices = geo.get_indices(geo.n_o+1);

            obj.f1 = f(obs_indices);
            f(obs_indices) = [];
            obj.f2= f(obs_indices);

            obj.f3= f(seg_indices);
            f(seg_indices) = [];
            obj.f4= f(seg_indices);
            
            obj.f5 = pd.gm*obj.f3;
            obj.f6 = pd.gm*obj.f4;
                
            
        end
    end
    
end