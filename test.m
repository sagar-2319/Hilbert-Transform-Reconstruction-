Afulll = Afull;
% A is your (2n+1)x(2n+1) matrix
n = (size(Afulll,1)-1)/2;
k = 500; % your desired value (must be <= n)

center = n + 1;
idxx = (center - k):(center + k);

A_restricted = Afulll(idxx, idxx);


number_of_modes_of_mu_one_side = 10000;

HT_matrix = A_restricted;

double_size = size(HT_matrix);
HT_matrix_size_one_side = (double_size(1)-1)/2;


%HT_matrix_size_one_side = m;
number_of_modes_of_L_one_side = 240;
number_of_modes_of_a_one_side = 700;
number_of_modes_of_h_one_side = 250;
resolution = 1;


%general rule for now - no_modes_h  <  no_modes_a ,  no_of_modes_L  <<
%HT_matrix_one_side, ideall half of the one sided row of HT, so if HT dim
%is 61x61, L should not have more than 15 modes


[mu_0,mu_p, mu_n] = mu_diff_ree(); 


[K_0,K_p,K_n] = K_diff_re(HT_matrix_size_one_side,HT_matrix);

[L_0, L_p, L_n,L_0_func,L_p_func,L_n_func] = L_diff_re(mu_0,mu_p,mu_n,K_0,K_p,K_n,number_of_modes_of_L_one_side, HT_matrix_size_one_side);
%[L_0, L_p, L_n,L_0_func,L_p_func,L_n_func] = L_diff_inf_re(mu_0,mu_p,mu_n,K_0,K_p,K_n,number_of_modes_of_L_one_side, HT_matrix_size_one_side,number_of_modes_of_mu_one_side);

[L_0_integral,L_p_integral,L_n_integral,E_0_func,E_p_func,E_n_func] = E_diff_re(number_of_modes_of_L_one_side,L_0,L_p,L_n);


[a,a_0,a_p,a_n,a_vec,theta,theta_prime,a_0_integrand_1,a_p_integrand_1,a_n_integrand_1] = a_coefficients_finufft_large_matlab(L_0_integral,L_p_integral,L_n_integral,L_0_func,L_p_func,L_n_func,E_0_func,E_p_func,E_n_func,number_of_modes_of_a_one_side);

[h,k_func,phi] = reconstruction_diff_re_large(a,number_of_modes_of_h_one_side,resolution,theta_prime);


