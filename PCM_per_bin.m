%% PCM IMPLEMENTATION

%% per bin!

% input: SH coefficients alpha, supported order N_sup,
% sharpness parameter eta

% output: S_hat, Gamma_sym, Gamma_norm, Gamma_phase, Gamma_S

% N_sup depends on analysis band (ex. 100 - 4000 Hz) and array
% configuration

% alpha(n) belongs to C^(2n+1)

% alpha = [alpha_00,         0,        0,      0, ..., 0,
%          alpha_1,-1,     alpha_1,0, alpha_11, 0, ..., 0,
%          ............
%          alpha_Nmax,-Nmax,      ...,   alpha_Nmax_Nmax ]

% rows: N+1
% columns: 2N+1

function [S_hat, G_sym, G_norm, G_phase, G_tot, G_tot_glob] = PCM_per_bin(alpha, N_sup, eta)
    S = 0;

    for j=1:N_sup+1
        n = j-1;
        alpha_norm(j) = norm(alpha(j,:));
        S(j) = alpha_norm(j)/sqrt(4*pi*(2*n+1));
    end
    S_hat_mod = min(S);
    
    % phase
    w0 = 1;
    alpha_00 = alpha(1,1);
    phi_0_cand = angle(alpha_00);
    delta_n0   = zeros(N_sup+1);
    w_n0        = zeros(N_sup+1);
    
    phi_n0_cand = zeros(N_sup+1);
    phi_nm_cand = zeros(N_sup+1,2*N_sup+1);
    eps = 1e-12;
    
    w_nm = zeros(N_sup+1, 2*N_sup+1);

    
    for j=1:N_sup+1
        n=j-1; % ordine reale
        phi_n0_cand(j) = angle(alpha(j,j)*(1i^(-n)));
        n1 = wrapToPi(phi_n0_cand(j) - phi_0_cand);
        n2 = wrapToPi(phi_n0_cand(j) + pi - phi_0_cand);
        if abs(n1) > abs(n2)
            phi_n0_cand(j) = phi_n0_cand(j) + pi;
        end
        delta_n0(j) = min(abs(n1),abs(n2));
    
        %reminder: coeff. (n,-m) is (j,j-m) and (n,m) is (j,j+m) ??in teoria??
        for k = 1:j % qua non serve che parta da 0 perché gli n0 li ha già calcolati
            m = k-j;
            rho_nm = conj(alpha(j,j-m))/((-1)^(n+m)*alpha(j,(j+m)) + eps);
            delta_nm = abs(log(abs(rho_nm) + eps));
            b = 2*abs(alpha(j,j+m))*abs(alpha(j,j-m))/((abs(alpha(j,j+m)))^2 + (abs(alpha(j,j-m)))^2 + eps);
            w_nm(j,k) = b/(1+delta_nm);
            phi_nm_cand(j,k) = -(1/2)*angle(rho_nm);
            n1 = wrapToPi(phi_nm_cand(j,k) - phi_0_cand);
            n2 = wrapToPi(phi_nm_cand(j,k) + pi - phi_0_cand);
            if abs(n1) > abs(n2)
                phi_nm_cand(j,k) = phi_nm_cand(j,k) + pi;
            end
        end
    end
    
    for j=1:N_sup+1
        w_n0(j) = (abs(alpha(j,j))/(sum(abs(diag(alpha))) + eps))*(cos(delta_n0(j)))^2;
    end
    
    sum_n0 = sum(w_n0(:) .* exp(1i * phi_n0_cand(:)));
    sum_nm = 0;
    for j = 2:N_sup+1
        n = j-1;
        for k = 1:j
            m = k-j;
            sum_nm = sum_nm + w_nm(j,k) * exp(1i * phi_nm_cand(j,k));
        end
    end
    
    phi_est = angle(1 * exp(1i * phi_0_cand) + sum_n0 + sum_nm);
    S_hat = S_hat_mod * exp(1i * phi_est);
    
    d_hat = alpha./(S_hat + eps);
    
    
    % così d_hat dovrebbe avere la stessa struttura di alpha, d_n sarebbero le
    % righe di d, ognuna che parte da n,-n e arriva a n,n + il resto della riga
    % con zero padding
    
    %% residuals
    % e_sym, e_norm, e_phase
    
    e_sym = zeros(1,N_sup+1);
    e_norm = zeros(1,N_sup+1);
    e_phase = zeros(1,N_sup+1);
    
    G_sym = zeros(1,N_sup+1);
    G_norm = zeros(1,N_sup+1);
    G_phase = zeros(1,N_sup+1);
    G_tot = zeros(1,N_sup+1);
    
    for j=1:N_sup+1
        n = j-1;
        e_norm(j) = abs(alpha_norm(j) / ((S_hat_mod + eps)*sqrt(4*pi*(2*n+1))) - 1);
        G_norm(j) = exp(-eta*e_norm(j));
        n1 = wrapToPi(angle(d_hat(j,j)) - n*pi/2);
        n2 = wrapToPi(angle(d_hat(j,j)) - pi - n*pi/2);
        e_phase(j) = (2/pi) * min(abs(n1), abs(n2));
        G_phase(j) = exp(-eta*e_phase(j));
        if n > 0
            sum_num = 0;
            sum_den = 0;
            for k = 1:j
                m = k-j;
                sum_num = sum_num + abs(d_hat(j, j-m) - (-1)^(n+m) * conj(d_hat(j, j+m)));
                sum_den = sum_den + abs(d_hat(j, j-m)) + abs(d_hat(j, j+m));
            end
            e_sym(j) = sum_num / (sum_den + eps);
            G_sym(j) = exp(-eta*e_sym(j));
        else
            e_sym(j) = 0;
        end
    
        G_tot(j) = G_sym(j).*G_norm(j).*G_phase(j);
    end
    
    if N_sup > 0
        G_tot_glob = mean(G_tot(2:N_sup+1));
    else
        G_tot_glob = G_tot(1);
    end
end
