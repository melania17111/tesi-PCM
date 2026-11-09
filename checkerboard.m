%% checkerboard test
% la struttura è circa quella del main

N_max = 3;
M = (N_max+1)^2;
SNR = 40;
eta = 2;
r = 0.04;
fs = 16000;
c = 343;

theta3 = pi/2;
phi3 = pi/4;
theta2 = pi/2;
phi2 = -pi/3;
theta1 = pi/4;
phi1 = pi;
T = linspace(0, 180, 181);
F = linspace(100, 5000, 200);

%% alpha

alpha = zeros(length(F), length(T), M);

for i=1:length(F)
    S_mag = 1;
    S1_phase = 2*pi*rand;
    S2_phase = 2*pi*rand;
    S3_phase = 2*pi*rand;
    S1 = S_mag*exp(1i*S1_phase);
    S2 = S_mag*exp(1i*S2_phase);
    S3 = S_mag*exp(1i*S3_phase);
    for j=1:length(T)
        if T(j)<=60
            alpha(i,j,:) = reshape(S1.*compute_d_plane(N_max,theta1,phi1), 1, 1, []);
        elseif T(j)<=120
            alpha(i,j,:) = reshape((S1.*compute_d_plane(N_max, theta1, phi1) + S2.*compute_d_plane(N_max, theta2, phi2))/sqrt(2), 1, 1, []);
        else
            alpha(i,j,:) = reshape((S1.*compute_d_plane(N_max, theta1, phi1) + S2.*compute_d_plane(N_max, theta2, phi2) + S3.*compute_d_plane(N_max, theta3, phi3))/sqrt(3), 1, 1, []);
        end
    end
end

% noise
sigP = mean(abs(alpha(:)).^2);
SNRlin = 10^(SNR/10);
noiseP = sigP / SNRlin;
noise_std = sqrt(noiseP/2);
noise = noise_std*(randn(size(alpha)) + 1i*randn(size(alpha)));
alpha = alpha + noise;

%% PCM

% preallocate
G_sym = zeros(length(F),length(T), N_max+1);
G_norm = zeros(length(F),length(T), N_max+1);
G_phase = zeros(length(F),length(T), N_max+1);
G_tot = zeros(length(F),length(T), N_max+1);
G_tot_glob = zeros(length(F),length(T));

% per TF bin now
for i=1:length(F)
    N_sup = min(N_max, floor(2*pi*r*F(i)/c));
    fprintf("N_sup nel main")
    disp(N_sup)
        for j=1:length(T)
            alpha_per_bin = reshape(alpha(i,j,:), 1, []);
            alpha_mat_per_bin = alpha_vector_to_matrix(alpha_per_bin, N_sup);
            fprintf("tempo")
            disp(T(j))
            [~, G1, G2, G3, G4, G5] = PCM_per_bin(alpha_mat_per_bin, N_sup, eta);
    
            G_sym(i,j,1:numel(G1))  = reshape(G1, 1, 1, []);
            G_norm(i,j,1:numel(G2)) = reshape(G2, 1, 1, []);
            G_phase(i,j,1:numel(G3))= reshape(G3, 1, 1, []);
            G_tot(i,j,1:numel(G4))  = reshape(G4, 1, 1, []);
            G_tot_glob(i,j) = G5;
        end
end

%%

figure(2)
for i=1:N_max+1
    subplot(2,2,i)
    imagesc(20*log10(G_tot(:,:,i)));
    axis xy;
    clim([-20 0])
    xlabel('T'); ylabel('F');
    colorbar;
end
sgtitle('G_{tot} for 1, 2 and 3 plane waves');

% le differenze tra la prima sezione (0<t<60) e le altre due ci sono ma non sono
% così nette come mi aspetterei

%%
figure(3)
imagesc(20*log10(G_tot_glob));
axis xy;
clim([-20 0])
xlabel('T'); ylabel('F');
colorbar; title('G_{tot glob}');