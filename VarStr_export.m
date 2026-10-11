% Plain-text export of VarStr.mlx (open the .mlx in MATLAB for the outputs)
function [Y,X] = VarStr(y,c,p)
% Build VAR regression matrices
% y: T x n
% c: 1 constant, 0 no constant
% p: lag order

[T,n] = size(y);

Y = y(p+1:end,:);                 % (T-p) x n
Teff = size(Y,1);

X = [];
if c==1
    X = ones(Teff,1);
end

for j = 1:p
    X = [X, y(p+1-j:end-j,:)];    % (T-p) x n each lag block
end
end
