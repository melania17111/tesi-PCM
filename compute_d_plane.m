function [d_plane] = compute_d_plane(N_max, theta, phi)

d_plane = zeros(1, (N_max+1)^2);
idx = 0;

for n = 0:N_max
    for m = -n:n
        Pnm_all = legendre(n, cos(theta));
        if n~=0
            Pnm = squeeze(Pnm_all(abs(m)+1));
        else
            Pnm = Pnm_all(1);
        end
        Nnm = sqrt((2*n+1)*factorial(n-abs(m))/(4*pi*factorial(n+abs(m))));
        Ynm = Nnm*Pnm*exp(1i*abs(m)*phi);
        if m<0
            Ynm = (-1)^m*conj(Ynm);
        end
        dnm = 4*pi*(1i^n)*conj(Ynm);
        idx = idx + 1;
        d_plane(idx) = dnm;
    end
end