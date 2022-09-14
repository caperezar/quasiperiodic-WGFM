function e = POU(t,t0,t1)
% Function used to construct the partition of the 
% unity. Based on 
% "Surface scattering in three dimensions: an accelerated
% high-oder solver" O. P. Bruno and L. A. Kunyansky
e = zeros(size(t));
for n=1:numel(t)
    if abs(t(n))<=t0
        e(n) = 1;
    elseif t0<abs(t(n)) && abs(t(n))<t1
        x = (abs(t(n))-t0)/(t1-t0); 
        e(n) = exp(2*exp(-1/x)/(x-1));
    elseif abs(t(n))>=t1
        e(n) = 0;
    end
end
   