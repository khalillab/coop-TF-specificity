clear all, clc

%% Initialize

TF = 10^0;                          % TF Concentration (µM)
Kd_r = logspace(-2,2,100)';         % DNA binding affinity (µM)
c_r = linspace(0,20,100)';          % Cooperativity (KbT)

% Max Txn rate
r0 = 1;                       

% # of binding sites
N1 = 1;      % Single
N4 = 4;      % Cooperative


%% Selectivity Parameter Space

% Single binding site dose response
txn_N1 = txn_fun( TF, Kd_r, N1, r0 );

% Initialize output
specificity = zeros(100,100);

% Calculate transcription for n=4 binding sites across parameters
for i=1:length(c_r)
    txn_N4 = txn_coops_fun( TF, Kd_r, N4, c_r(i), r0);
    
    % Specificity score
    specificity(i,:) = txn_N4-txn_N1;
end


% Define color map for plotting
mymap = [255,247,243;
253,224,221;
252,197,192;
250,159,181;
247,104,161;
221,52,151;
174,1,126;
122,1,119;
73,0,106]/255;


% Plot parameter space
figure
    contourf(log10(Kd_r),c_r,specificity,8)
    colormap(mymap)
    pbaspect([1 1 1])
    set(gca,'FontSize',18)
    colorbar
    xlabel('affinity (K_T_F)')
    ylabel('cooperativity (c)')

    
    
%% Dose response Parameters

% TF Range
TF_r = logspace(-3,0,100)';

% Affinity terms (Ktf)
Kd_low  = 10^-1;        % high affinity
Kd_med  = 10^0;         % medium affinity
Kd_high = 10^1;         % low affinity

% Cooperativity terms (c)
c_low = 2.5;            % low coop
c_med = 10;             % medium coop
c_high = 20;            % high coop


%% Plot Dose Responses
figure

    % Low Ktf, Low c
    subplot(3,1,1)
    semilogx(TF_r,txn_fun( TF_r, Kd_low, N1, r0 ),'-'); hold on
    semilogx(TF_r,txn_coops_fun( TF_r, Kd_low, N4, c_low, r0 ),'-'); hold on
    ylim([0 1.1])
    xlim([min(TF_r) max(TF_r)])
    pbaspect([1 1 1])
    
    % Medium Ktf, Medium c
    subplot(3,1,2)
    semilogx(TF_r,txn_fun( TF_r, Kd_med, N1, r0 ),'-'); hold on
    semilogx(TF_r,txn_coops_fun( TF_r, Kd_med, N4, c_med, r0 ),'-'); hold on
    ylim([0 1.1])
    xlim([min(TF_r) max(TF_r)])
    pbaspect([1 1 1])
    
    % High Ktf, High c
    subplot(3,1,3)
    semilogx(TF_r,txn_fun( TF_r, Kd_high, N1, r0 ),'-'); hold on
    semilogx(TF_r,txn_coops_fun( TF_r, Kd_high, N4, c_high, r0 ),'-'); hold on
    ylim([0 1.1])
    xlim([min(TF_r) max(TF_r)])
    pbaspect([1 1 1])
    xlabel('[TF]')
    ylabel('transcriptional output')

    legend('txn (n=1)','txn (n=4)')