% Plain-text export of SignZeroRestrictBoot.mlx (open the .mlx in MATLAB for the outputs)
function [med, band, acc] = SignZeroRestrictBoot(y,p,H,c,MaxBoot,cl,k0,lrestr,rest_z)
% Identificazione Sign-Zero (HM2 style)
% k0: n. variabili slow (zero impact)

[T,n] = size(y);
m = n - k0; % dimensione blocco fast

if length(rest_z) ~= m
    error("rest_z deve essere (n-k0) x 1");
end

% VAR iniziale e matrici base
[Y,X] = VarStr(y,c,p);
Bols  = (X'*X)\(X'*Y);
u     = Y - X*Bols;
B     = companion(reshape(Bols, (c+n*p)*n, 1),p,n,c);

store = nan(n,H,MaxBoot);
acc   = 0;
uT    = u'; 

for b = 1:MaxBoot
    % --- Bootstrap serie ---
    Teff = size(Y,1);
    YB = zeros(n*p, Teff+1);
    YB(:,1) = X(1,2:end)'; 
    for t = 2:Teff+1
        uu = uT(:, randi(size(uT,2))); 
        YB(:,t) = [Bols(1,:)'; zeros(n*p-n,1)] + B*YB(:,t-1) + [uu; zeros(n*p-n,1)];
    end
    yb = [y(1:p,:); YB(1:n,2:end)'];

    % --- Stima su campione boot ---
    [Yb,Xb] = VarStr(yb,c,p);
    Bb_ols  = (Xb'*Xb)\(Xb'*Yb);
    Sb      = chol(cov(Yb - Xb*Bb_ols), 'lower');
    Bb_comp = companion(reshape(Bb_ols, (c+n*p)*n, 1),p,n,c);

    % Pre-calcolo risposte Cholesky
    phi_s = zeros(n,n,H);
    for h=1:H
        tmp = Bb_comp^(h-1);
        phi_s(:,:,h) = tmp(1:n,1:n) * Sb;
    end

    % --- Prova rotazioni ---
    for tr=1:20000
        Q = RotationMatrix(m, randn(m*(m-1)/2, 1));
        
        % Vettore rotazione: identità sopra, Q sotto
        q_col = [zeros(k0,1); Q(:,1)]; 
        
        itmp = zeros(n,H);
        for h=1:H
            itmp(:,h) = phi_s(:,:,h) * q_col;
        end

        % Verifica segni sul blocco fast
        if all(sign(itmp(k0+1:end, lrestr)) == rest_z)
            acc = acc + 1;
            store(:,:,b) = itmp;

            % Stampa matrice H per controllo (solo al primo colpo)
            if acc == 1
                fprintf('\n--- Matrice H (Identity + Rotation) ---\n');
                H_full = eye(n);
                H_full(k0+1:end, k0+1:end) = Q;
                disp(round(H_full, 4));
                fprintf('---------------------------------------\n');
            end
            break;
        end
    end
end

% Statistiche finali
med  = median(store,3,'omitnan');
band = prctile(store,cl,3); 

end
