% Plain-text export of Project3.mlx (open the .mlx in MATLAB for the outputs)
close all
clear all
clc

data = readmatrix("excelprogetto.xlsx","Sheet","Monthly","Range","B2:K770");

% y = IP, Infl, FedFunds, Spread, SP500
y = [100*diff(log(data(:,3))) , 100*diff(log(data(:,2))) , data(2:end,1) , data(2:end,6) , 100*diff(log(data(:,8)))];

high = data(2:end,10);

vNames = {'IP growth','Inflation','FedFunds','BAA-AAA','SP500 ret'};
shkIdx = 3; % shock MP su FedFunds (per Cholesky)

% Split dei campioni per regime
y_h = y(high==1,:);
y_l = y(high==0,:);

%% optLag (AIC)
max_p = 12;
AIC   = nan(max_p,1);
n = size(y,2);

for p = 1:max_p
    [~,res] = VARTopicsOLS(y,p);
    Teff  = size(res,1);
    Sigma = (res'*res)/Teff;
    k     = n^2*p + n;
    AIC(p)= log(det(Sigma)) + 2*k/Teff;
end

[~,popt] = min(AIC);
fprintf("Lag scelto = %d\n", popt);

%% IRF e Bootstrap (Cholesky)
H = 36;
c = 1;
MaxBoot = 1000;
cl = [16 84];
ind = [];
opt = 0;

[cirfL,CholBootL] = CholeskyBoot(y_l,popt,H,c,MaxBoot,cl,ind,opt);
[cirfH,CholBootH] = CholeskyBoot(y_h,popt,H,c,MaxBoot,cl,ind,opt);

% Risposta a shock MP (Cholesky)
irfL = squeeze(cirfL(:,shkIdx,:));
irfH = squeeze(cirfH(:,shkIdx,:));

% Percentili bande (Cholesky)
tmpL = prctile(CholBootL(:,shkIdx,:,:),cl,4);
tmpH = prctile(CholBootH(:,shkIdx,:,:),cl,4);

bandL = permute(squeeze(tmpL),[1 3 2]);
bandH = permute(squeeze(tmpH),[1 3 2]);

%% plot cholesky
t = 1:H;
figure('Name','IRF Cholesky: Low vs High','Color','w');
for i=1:n
    subplot(n,1,i); hold on;
    plot(t, squeeze(bandL(i,1,:)), 'b:', 'LineWidth', 0.9);
    plot(t, squeeze(bandL(i,2,:)), 'b:', 'LineWidth', 0.9);
    plot(t, squeeze(bandH(i,1,:)), 'r:', 'LineWidth', 0.9);
    plot(t, squeeze(bandH(i,2,:)), 'r:', 'LineWidth', 0.9);
    plot(t, irfL(i,:), 'b', 'LineWidth', 1.3);
    plot(t, irfH(i,:), 'r', 'LineWidth', 1.3);

    yline(0,'k-');
    grid on; xlim([1 H]);
    title(vNames{i},'Interpreter','none');
    if i==1
        legend('Low','','High','','Low IRF','High IRF','Location','best');
    end
end

%% SIGN–ZERO RESTRICTIONS


kZero   = 2;          % IP, Inflation are slowmoving
lagrest = 1;          % impact horizon

% Signs for FFR, Spread, SP500
rest_z = [+1; +1; -1];  

% LOW
[medL_z, bandL_z, accL] = SignZeroRestrictBoot(y_l,popt,H,c,MaxBoot,cl,kZero,lagrest,rest_z);
% HIGH
[medH_z, bandH_z, accH] = SignZeroRestrictBoot(y_h,popt,H,c,MaxBoot,cl,kZero,lagrest,rest_z);

fprintf("Zero-sign accepted draws: LOW=%d / %d, HIGH=%d / %d\n", accL, MaxBoot, accH, MaxBoot);

figure('Name','Sign-Zero Restrictions: Low vs High','Color','w');
for i=1:n
    subplot(n,1,i); hold on;
    plot(t, squeeze(bandL_z(i,:,1)), 'b:', 'LineWidth', 0.9);
    plot(t, squeeze(bandL_z(i,:,2)), 'b:', 'LineWidth', 0.9);
    plot(t, squeeze(bandH_z(i,:,1)), 'r:', 'LineWidth', 0.9);
    plot(t, squeeze(bandH_z(i,:,2)), 'r:', 'LineWidth', 0.9);
    plot(t, medL_z(i,:), 'b', 'LineWidth', 1.3);
    plot(t, medH_z(i,:), 'r', 'LineWidth', 1.3);

    yline(0,'k-'); grid on; xlim([1 H]);
    title(vNames{i},'Interpreter','none');
    if i==1
        legend('Low','','High','','Low med','High med','Location','best');
    end
end

%% FEVD (Cholesky)
hor = [12 24 36];
shkNames = matlab.lang.makeValidName(vNames);

[~, vdkL] = TopicsVD2(cirfL, cirfL, hor);
[~, vdkH] = TopicsVD2(cirfH, cirfH, hor);

for k = 1:numel(hor)
    ML_k = squeeze(vdkL(:,k,:));
    MH_k = squeeze(vdkH(:,k,:));
    MD_k = MH_k - ML_k;

    TL = array2table(round(ML_k,2), 'RowNames', vNames, 'VariableNames', shkNames);
    TH = array2table(round(MH_k,2), 'RowNames', vNames, 'VariableNames', shkNames);
    TD = array2table(round(MD_k,2), 'RowNames', vNames, 'VariableNames', shkNames);

    fprintf("\nFEVD @ %d mesi\n", hor(k));
    disp("LOW");  disp(TL);
    disp("HIGH"); disp(TH);
    disp("DIFF"); disp(TD);
end
