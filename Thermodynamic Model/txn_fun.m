function [ txn ] = txn_fun( TF, Ktf, N, r0 )

% Fraction bound
fb = [];

% Txn rates
ratesTF = r0*[1:N]/N;

%% TF-DNA States
for n=1:N
    mult = nchoosek(N,n);            % Multiplicity
    fb  = [fb mult*(TF./Ktf).^n];
end

%% Calculate Mean Transcription
rates = ratesTF';
txn = (fb*rates)./(1 + fb*ones(length(rates),1));


end

