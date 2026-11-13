function alpha_complex = alpha_real_to_complex(alpha, N_max)
    alpha_complex = zeros(size(alpha));
    for j=1:N_max+1
        n = j-1;
        for m = -n:n
            if m==0
                Y_nm = alpha(j,j+m);
            elseif m > 0
                Y_nm = (1/ sqrt(2)) .* (alpha(j,j+m) - 1i*alpha(j,j-m));
            else
                Y_nm = ((-1)^m / sqrt(2)) .* (alpha(j,j+abs(m)) + 1i*alpha(j,j-(abs(m))));
            end
            alpha_complex(j,j+m) = (1i^n) .* Y_nm;
        end
    end