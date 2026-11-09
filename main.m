% to be given: input signal, N_sup, eta

%% usare angle
% dataset di 50 stanze, snr sempre 40, cambiano dimensioni, assorbimento e
% posizione di microfono e sorgente

%%
filename = 'dataset_stanze_ottimizzato.h5';
room_idx = 5; % indice della stanza da analizzare
group_path = sprintf('/room_%d', room_idx);

raw_rirs = h5read(filename, [group_path '/rirs']);

hoa_rirs = raw_rirs(:,:,1);

positions   = h5read(filename, [group_path '/mic_positions']);
positions = positions(:,1);
source_pos = h5read(filename, [group_path '/source_pos']);
room_dim   = h5read(filename, [group_path '/room_dim']);
absorption  = h5readatt(filename, group_path, 'absorption');

fprintf('Stanza %d caricata correttamente: RIR di dimensione [%d x %d]\n', ...
    room_idx, size(hoa_rirs, 1), size(hoa_rirs, 2));

fs = 16000;
c = 343;

%%
%SNR = snr_db;
N_max = 3;
eta = 1;
r = 0.04;
fs = double(fs);
[source, fs_s] = audioread('LJ001-0001.wav'); %speech
source = source(1:50000);
disp(length(source))
[p, q] = rat(fs/fs_s);
source = resample(source, p, q);
source = source/(max(source));
disp(length(source))
[alpha, F, T] = compute_alpha(hoa_rirs, source, fs, N_max);
fprintf("F è lunga")
disp(length(F))
M = (N_max+1)^2;

%%
nfft = 512;
audio = istft(alpha(:,:,1), fs, 'Window', hann(nfft), ...
    'OverlapLength', nfft/2, 'FrequencyRange', 'onesided');

plot(audio)
%% e poi
% preallocate
G_sym = zeros(length(F),length(T), N_max+1);
G_norm = zeros(length(F),length(T), N_max+1);
G_phase = zeros(length(F),length(T), N_max+1);
G_tot = zeros(length(F),length(T), N_max+1);
G_tot_glob = zeros(length(F),length(T));
alpha_mat_per_bin = zeros(N_max+1, 2*N_max+1);
alpha_per_bin = zeros(1, (N_max+1)^2);
% F x T x M

% bisogna salvare alpha per ogni frequenza x fare la ricostruzione e il check
alpha_total = zeros(length(F), (N_max+1)^2);

% per TF bin now
for i=1:length(F)
    N_sup = min(N_max, floor(2*pi*r*F(i)/c));
    fprintf("check: N_sup nel main")
    disp(N_sup)
    alpha_avg = 0;
    for j=1:length(T)
        alpha_per_bin = reshape(alpha(i,j,:), 1, []);
        alpha_mat_per_bin = alpha_vector_to_matrix(alpha_per_bin, N_sup);
        
        [~, G1, G2, G3, G4, G5] = PCM_per_bin(alpha_mat_per_bin, N_sup, eta);

        G_sym(i,j,1:numel(G1))  = reshape(G1, 1, 1, []);
        G_norm(i,j,1:numel(G2)) = reshape(G2, 1, 1, []);
        G_phase(i,j,1:numel(G3))= reshape(G3, 1, 1, []);
        G_tot(i,j,1:numel(G4))  = reshape(G4, 1, 1, []);
        G_tot_glob(i,j) = G5;
        alpha_avg = alpha_avg + alpha_per_bin;
    end
    alpha_total(i, :) = alpha_avg/(length(T));
end

%% check

RES = [50 50];
elev = linspace(-pi/2, pi/2,RES(1)); % inclination
phi   = linspace(0,2*pi,RES(2)); % azimuth
Ynm = zeros(M, length(elev), length(phi));

%%
idx = 1;
theta = -elev + pi/2;
for n = 0:N_max
    P_all = legendre(n, cos(theta));
    for m = -n:n
        m_abs = abs(m);
        Plm = reshape(P_all(m_abs+1, :), [1, length(theta)]); % 1 x ntheta
        Nnm = sqrt((2*n+1)/(4*pi) * factorial(n-m_abs)/factorial(n+m_abs));

        [Phi, Theta] = meshgrid(phi, theta);
        Plm2 = repmat(Plm(:), 1, length(phi)); % ntheta x nphi
        Ypos = Nnm * Plm2 .* exp(1i * m_abs * Phi); % ntheta x nphi

        if m < 0
            Ycur = (-1)^m_abs * conj(Ypos);   % Y_l^{-m}
        else
            Ycur = Ypos;                      % Y_l^{m}
        end

        Ynm(idx,:,:) = Ycur;
        idx = idx + 1;
    end
end

%% pseudospectrum

N_coeffs = (N_max + 1)^2;
N_grid = length(theta) * length(phi);
N_TF = length(F) * length(T);

alpha_flat = reshape(alpha, N_TF, N_coeffs);
Y_flat = reshape(Ynm, N_coeffs, N_grid);

% calcolo di p_ft per tutti i bin TF e tutte le direzioni contemporaneamente
P_ft_all = alpha_flat * conj(Y_flat);

E_tot_flat = mean(abs(P_ft_all).^2, 1);

E_tot = reshape(E_tot_flat, length(theta), length(phi));

%% plot canali e balloon totale

figure(7)
E_tot_norm = E_tot / max(E_tot(:));
plot3dBalloon(E_tot_norm, phi*180/pi, elev*180/pi, positions, 'Xlim', [0 room_dim(1)], 'Ylim', [0 room_dim(2)], 'Zlim', [0 room_dim(3)])
hold on
scatter3(source_pos(1), source_pos(2), source_pos(3), '*')

%% plot G

figure(2)
for i=1:N_max+1
    subplot(2,2,i)
    imagesc(T(1:end), F(1:end), 20*log10(G_norm(:,:,i)));
    axis xy;
    clim([-10, 0])
    xlabel('T'); ylabel('F');
    colorbar; title('G_{tot}');
end

% g sym 1 è al contrario
% g phase non si capisce niente
% g norm confuso

%% plot G glob e confronto con spettrogramma sorgente
figure(8)
x1 = subplot(2,1,1);
imagesc(T(1:end), F(1:end), 20*log10(G_tot_glob));
axis xy;
clim([-15 0])
xlabel('T'); ylabel('F');
%colorbar;
title('G_{tot glob}');

nfft = 512;
x2 = subplot(2,1,2);
spectrogram(source, hann(nfft), nfft/2, nfft, fs, 'yaxis');
[S, F_s, T_s] = spectrogram(source, hann(nfft), nfft/2, nfft, fs);
imagesc(T_s, F_s, 20*log10(abs(S)));
axis xy;
linkaxes([x1, x2], 'x');
title('source stft for comparison')


%% plot hoa

figure(4)
%for i=1:M
 %   subplot(4,4,i)
    plot(squeeze(raw_rirs(:,:,1)))
%end
sgtitle('RIRS for each order and grade (abs)')

%% source
figure(9)
plot(source)

%%
mic_pos = positions;
