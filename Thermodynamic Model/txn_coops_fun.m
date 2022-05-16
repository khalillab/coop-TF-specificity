function [ txn ] = txn_coops_fun( TF, Ktf, N, c, r0 )

% Fraction bound
fb = [];

% Txn rates
ratesTF = r0*[1:N]/N;
ratesC = r0;

%% TF-DNA States
for n=1:N
    deg = nchoosek(N,n);        % Degeneracy
    fb  = [fb deg*(TF./Ktf).^n];
end

rates = [ratesTF'];


%% TF-DNA-Clamp State
if c~=0
    rates = [ratesTF'; ratesC'];
    fb  = [fb ((TF./Ktf).^n)/exp(-c)];
end


%% Calculate Mean Transcription
txn = (fb*rates)./(1 + fb*ones(length(rates),1));


end

