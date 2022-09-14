

function [phi_1,psi_1,phi_2,psi_2] = traces(mu,geo)

obs_indices = geo.get_indices(1);
seg_indices = geo.get_indices(2);

phi_1 = mu(obs_indices);
mu(obs_indices) = [];
psi_1 = mu(obs_indices);
  

phi_2 = mu(seg_indices);
mu(seg_indices) = [];
psi_2 = mu(seg_indices);