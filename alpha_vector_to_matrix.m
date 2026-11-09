function [alpha_mat] = alpha_vector_to_matrix(alpha, N)
alpha_mat = zeros(N+1, 2*N+1);
% qua non passa alpha intera ma solo al bin ij quindi alpha(i,j,:)
counter = 1;
for n = 1:N+1
    alpha_mat(n,1:numel(alpha(counter:(counter+2*(n-1))))) = alpha(counter:(counter+2*(n-1)));
    counter = counter + 2*(n-1) + 1;
end
