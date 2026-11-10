%% extraction of alpha_nm
% importo le h_nm dallo script python
% inizio con un microfono solo nel far field
function [alpha, F, T] = compute_alpha(h_single, source, fs, N_max)
M = (N_max+1)^2;
nfft = 512;
source = source(:);

num_samples = length(source) + size(h_single,1) - 1;
x_hoa = zeros(num_samples, M);

for j = 1:M
    x_hoa(:, j) = conv(source, h_single(:, j));
end

[S_full, F_full, T] = stft(x_hoa(:,1), fs, 'Window', hann(nfft), ...
    'OverlapLength', nfft/2, 'FrequencyRange', 'centered');

% freq positive
pos_idx = (F_full >= 0);
F = F_full(pos_idx);
disp(size(F_full))
disp(size(F))
alpha = zeros(length(F), length(T), M);
alpha(:,:,1) = S_full(pos_idx, :);

for j = 2:M
    [S_temp, ~, ~] = stft(x_hoa(:,j), fs, 'Window', hann(nfft), ...
        'OverlapLength', nfft/2, 'FrequencyRange', 'centered');
    alpha(:,:,j) = S_temp(pos_idx, :);
end
end
% da qua ottengo alpha FxTxM quindi poi mi serve una funzione per
% convertirla nella matrice che ho definito all'inizio di PCM.m