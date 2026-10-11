function R = RotationMatrix(q, theta) % q è la dimensione della rotation matrix quindi deve essere uguale a n
% tetha è il vettore di parametri di dimensione n*(n-1)/2
R = eye(q);
M = R;
k = 1;
for j = 1 : q - 1
    for i = j + 1 : q
        M(j , j) = cos(theta(k));
        M(j , i) = sin(theta(k));
        M(i , j) = - sin(theta(k));
        M(i , i) = cos(theta(k));
        R = R*M;
        M = eye(q);
        k = k + 1;
    end
end