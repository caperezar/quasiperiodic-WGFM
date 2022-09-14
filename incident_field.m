function [u,du]= incident_field(x,problem_data)
% all the angles are measure with respect to the horizontal line

alpha = problem_data.angle;

k = problem_data.k(1);

d = [cos(alpha) sin(alpha)];

u = exp(1i*k*dot(x,d));

du(1) = 1i*k*cos(alpha)*exp(1i*k*dot(x,d));
du(2) = 1i*k*sin(alpha)*exp(1i*k*dot(x,d));



end