% Plain-text export of companion.mlx (open the .mlx in MATLAB for the outputs)
function B = companion(VecB,p,n,c)
k    = c + n*p;
Bols = reshape(VecB,k,n);         % (c+n*p) x n

A = Bols(c+1:end,:)';             % n x (n*p)

B = [A;
     eye(n*(p-1)) zeros(n*(p-1),n)];
end
