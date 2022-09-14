classdef physical_data
    
    properties 
        k
        nu
        gm
        polarization
        angle        
        L
        alpha
        beta
        
        alpha_n
        beta_n
        prop_modes 

        reg_alpha_n
        reg_beta_n
        reg_modes

        amp_sine_curve
        freq_sine_curve

    end
    
    methods
        
        function obj = physical_data(k,pol,angle,L)
            
            obj.k = k;
    
            obj.polarization = pol;
                        
            obj.angle = angle;
            
            obj.L = L;
            
            obj.gm = exp(1i*k(1)*cos(angle)*L);
            
            if strcmp(pol,'TE')
                
                obj.nu = 1;
                
            elseif strcmp(pol,'TM')
                
                obj.nu = k(1)^2/k(2)^2;
                
            end
            
            obj.alpha = k(1)*cos(angle);
            
            obj.beta =  k(1)*sin(angle);
            
            K = 2*pi/L;
            
            N_min = -ceil(abs(k(1)+obj.alpha)/K);
            N_max = ceil(abs(k(1)-obj.alpha)/K);
            N_range = (N_min:N_max);
            obj.alpha_n = obj.alpha + K*N_range;
            
            obj.beta_n = sqrt(k(1)^2-obj.alpha_n.^2);
            
            obj.beta_n = abs(real(obj.beta_n)) + 1i*abs(imag(obj.beta_n));
            
            indx = find(obj.alpha_n.^2<=k(1)^2);
            obj.alpha_n = obj.alpha_n(indx);
            obj.beta_n = obj.beta_n(indx);
            obj.prop_modes = N_range(indx);
            
            %% Modes to be considered in the integral equation formulation
%             obj.reg_alpha_n = obj.alpha_n;
%             obj.reg_beta_n = obj.beta_n;
%             obj.reg_modes = obj.prop_modes;


%             N_range = (2*N_min:2*N_max);
            N_range = (-100:100);

            obj.reg_alpha_n = obj.alpha + K*N_range;
            obj.reg_beta_n = sqrt(k(1)^2-obj.reg_alpha_n.^2);
                        
%             indx = union(indx,find(abs(obj.reg_beta_n)<=1*k(1)/2));
            indx = find(abs(obj.reg_beta_n)<=3/4*k(1));
% %             indx = intersect(find(abs(obj.reg_alpha_n)>=k(1)/2),indx);
% 
            obj.reg_alpha_n = obj.reg_alpha_n(indx);
            obj.reg_beta_n = obj.reg_beta_n(indx);
            obj.reg_modes = N_range(indx);
        end
                
    end
    
end
        
        
        

