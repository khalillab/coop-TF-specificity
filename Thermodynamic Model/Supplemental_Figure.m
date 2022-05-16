clear all, clc

%% Initialize

TF = 10^0;                          % TF Concentration (µM)
Kd_r = logspace(-2,2,100)';         % DNA binding affinity (µM)
c_r = linspace(0,20,100)';          % Cooperativity (KbT)

% Max Txn rate
r0 = 1;                       
    
    
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


    
%% Selectivity Parameter Space for different N (Txn(N) - Txn(n=1))

for N=3:5

    % Single binding site dose response
    txn_N1 = txn_fun( TF, Kd_r, 1, r0 );

    % Initialize output
    specificity = zeros(100,100);

    % Calculate transcription for n=4 binding sites across parameters
    for i=1:length(c_r)
        txn_coops = txn_coops_fun( TF, Kd_r, N, c_r(i), r0);
        
        % Specificity score
        specificity(i,:) = txn_coops-txn_N1;
    end

    % Plot parameter space
    figure
        contourf(log10(Kd_r),c_r,specificity,8)
        colormap(mymap)
        pbaspect([1 1 1])
        set(gca,'FontSize',18)
        colorbar
        title(['txn(n=' , int2str(N), ') - txn(n=1)'])
end


%% Alternative Specificity score - division (Txn(n=4) / Txn(n=1))


% Single binding site dose response
txn_N1 = txn_fun( TF, Kd_r, 1, r0 );

% Define outputs
specificity = zeros(100,100);

% Calculate transcription for n=4 binding sites across parameters
for i=1:length(c_r)
    txn_coops = txn_coops_fun( TF, Kd_r, 4, c_r(i), r0);
    specificity(i,:) = txn_coops./txn_N1;
end


figure
    contourf(log10(Kd_r),c_r,log10(specificity),8)
    colormap(mymap)
    pbaspect([1 1 1])
    set(gca,'FontSize',18)
    colorbar
    title(['txn(n=4) / txn(n=1)'])
