function b = right_hand_side(pd,geo)         
phi = rhs;
F  = [];
dF = [];
for n=1:geo.n_o
    phi = phi.cctor_local(pd,geo,n);
    F  = [F;phi.IF];
    dF = [dF;phi.dIF];
end
np_obs = geo.n_o*geo.obstacles(1).Np;

b = zeros(2*(geo.n_o*geo.obstacles(1).Np+geo.segments(1).Np),1);

b(1:2*np_obs) = [F;dF];


end



