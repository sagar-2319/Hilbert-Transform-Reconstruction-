function [K_0,K_p,K_n] = K_diff_re(HT_matrix_size_one_side,HT_matrix)

K_0 = zeros(1,2*HT_matrix_size_one_side+1);
K_p = zeros(1,2*HT_matrix_size_one_side+1);
K_n = zeros(1,2*HT_matrix_size_one_side+1);


for mode = -HT_matrix_size_one_side:HT_matrix_size_one_side

        K_0(1,idx(mode,HT_matrix_size_one_side)) = 1i.*(exp(1i.*2.*pi.*(-HT_matrix_size_one_side:HT_matrix_size_one_side)./3)-exp(-1i.*2.*pi.*(-HT_matrix_size_one_side:HT_matrix_size_one_side)./3))*HT_matrix(idx(-mode,HT_matrix_size_one_side),:).';
        K_p(1,idx(mode,HT_matrix_size_one_side)) = 1i.*(1-exp(-1i.*2.*pi.*(-HT_matrix_size_one_side:HT_matrix_size_one_side)./3))*HT_matrix(idx(-mode,HT_matrix_size_one_side),:).';
        tic;K_n(1,idx(mode,HT_matrix_size_one_side)) = 1i.*(1-exp(1i.*2.*pi.*(-HT_matrix_size_one_side:HT_matrix_size_one_side)./3))*HT_matrix(idx(-mode,HT_matrix_size_one_side),:).';toc;
end

function idx = idx(k, m)
    idx = k + m + 1;
end

end