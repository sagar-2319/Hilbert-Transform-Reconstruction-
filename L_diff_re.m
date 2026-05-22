
function [L_0, L_p, L_n,L_0_func,L_p_func,L_n_func] = L_diff_re(mu_0,mu_p,mu_n,K_0,K_p,K_n,number_of_modes_of_L_one_side, HT_matrix_size_one_side)

L_0 = zeros(1,2*number_of_modes_of_L_one_side+1);
L_p = zeros(1,2*number_of_modes_of_L_one_side+1);
L_n = zeros(1,2*number_of_modes_of_L_one_side+1);

for mode = -number_of_modes_of_L_one_side:number_of_modes_of_L_one_side


    L_0(1,idx(mode,number_of_modes_of_L_one_side)) = mu_0(mode - (-HT_matrix_size_one_side:HT_matrix_size_one_side)) * K_0(1,idx(-HT_matrix_size_one_side:HT_matrix_size_one_side,HT_matrix_size_one_side)).';
    L_p(1,idx(mode,number_of_modes_of_L_one_side)) = mu_p(mode - (-HT_matrix_size_one_side:HT_matrix_size_one_side)) * K_p(1,idx(-HT_matrix_size_one_side:HT_matrix_size_one_side,HT_matrix_size_one_side)).';
    tic;L_n(1,idx(mode,number_of_modes_of_L_one_side)) = mu_n(mode - (-HT_matrix_size_one_side:HT_matrix_size_one_side)) * K_n(1,idx(-HT_matrix_size_one_side:HT_matrix_size_one_side,HT_matrix_size_one_side)).';toc;
        

end    



inputss = -number_of_modes_of_L_one_side:number_of_modes_of_L_one_side;

L_0_func = @(s) sum(L_0.*(exp(1i.*s.*inputss))).*(s >= -pi/3 & s < pi/3);
L_p_func = @(s) sum(L_p.*(exp(1i.*s.*inputss))).*(s >= pi/3 & s <= pi);
L_n_func = @(s) sum(L_n.*(exp(1i.*s.*inputss))).*(s >= -pi & s < -pi/3);


function idx = idx(k, m)
    idx = k + m + 1;
end

end