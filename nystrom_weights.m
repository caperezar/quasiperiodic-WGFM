function R = nystrom_weights(prt)

N = prt.N;

t = pi/N*(0:2*N-1);

tpmintj = repmat(t,2*N,1)-repmat(t.',1,2*N);

temp_vec = zeros([1,1,N-1]);

temp_vec(1,1,:)=(1:N-1);

R=-2*pi/N*sum(cos(repmat(tpmintj,[1,1,N-1]).*repmat(temp_vec,[2*N,2*N,1]))./...
    repmat(temp_vec,[2*N,2*N,1]),3)-pi/(N^2)*cos(N*tpmintj);

Np = prt.Np;

R  = R(1:Np,1:Np);
end