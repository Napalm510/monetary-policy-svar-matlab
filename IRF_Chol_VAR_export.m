% Plain-text export of IRF_Chol_VAR.mlx (open the .mlx in MATLAB for the outputs)
function irf = IRF_Chol_VAR(para,sigm,p,H,shock_idx)

n = size(sigm,1);

% companion from para
B = [para(2:end,:)'; eye(n*(p-1)) zeros(n*(p-1),n)];
P = chol(sigm,'lower'); 

irf = zeros(n,H);
e = zeros(n,1); e(shock_idx)=1;

for h = 1:H
    Phi = B^(h-1);                         
    G = Phi(1:n,1:n) * P;                  
    irf(:,h) = G * e;                    
end
end
