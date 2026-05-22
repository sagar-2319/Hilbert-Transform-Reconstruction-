function [a,a_0,a_p,a_n,a_vec,theta,theta_prime,...
          a_0_integrand_1_n,a_p_integrand_1_n,a_n_integrand_1_n] = ...
          a_coefficients_finufft_large_matlab( ...
          L_0_integral, L_p_integral , L_n_integral , ...
          L_0_func,L_p_func,L_n_func, ...
          E_0_func,E_p_func,E_n_func, ...
          number_of_modes_of_a_one_side)


%--------------------------------------------

a_0_integrand_1_n = @(s) atan(sqrt(3).*(1-E_0_func(s))./(1+E_0_func(s)));
a_p_integrand_1_n = @(s) atan(sqrt(3).*(1-E_p_func(s))./(1+E_p_func(s)));
a_n_integrand_1_n = @(s) atan(sqrt(3).*(1-E_n_func(s))./(1+E_n_func(s)));

a_0_integrand_2_n = @(s) (L_0_func(s).*E_0_func(s))./(1-E_0_func(s)+(E_0_func(s)).^2);
a_p_integrand_2_n = @(s) (L_p_func(s).*E_p_func(s))./(1-E_p_func(s)+(E_p_func(s)).^2);
a_n_integrand_2_n = @(s) (L_n_func(s).*E_n_func(s))./(1-E_n_func(s)+(E_n_func(s)).^2);

%-------------------------------------------

a_0 = zeros(1,2*number_of_modes_of_a_one_side+1);
a_p = zeros(1,2*number_of_modes_of_a_one_side+1);
a_n = zeros(1,2*number_of_modes_of_a_one_side+1);

%------------------------------------------

theta = @(s) ...
      2.*a_0_integrand_1_n(s).*(s >= -pi/3 & s <  pi/3) ...
    + (2.*a_p_integrand_1_n(s)+2*pi/3).*(s >=  pi/3 & s <= pi) ...
    + (-2.*a_n_integrand_1_n(s)-2*pi/3).*(s >= -pi & s < -pi/3);

theta = @(s) arrayfun(theta,s);

theta_prime = @(s) sqrt(3)./2.*( ...
     -a_0_integrand_2_n(s).*(s >= -pi/3 & s <  pi/3) ...
     +a_p_integrand_2_n(s).*(s >=  pi/3 & s <= pi) ...
     -a_n_integrand_2_n(s).*(s >= -pi & s < -pi/3) );

theta_prime = @(s) arrayfun(theta_prime,s);


%-------------------------------------------

tic;

M = number_of_modes_of_a_one_side;
N = 2^13; %2^18

s = linspace(-pi, pi, N+1);
s(end) = [];

theta_vals   = theta(s);
theta_p_vals = theta_prime(s);

t = real(theta_vals(:));
t = mod(t + pi, 2*pi) - pi;

w = complex((theta_p_vals(:)).^2);

t = double(t(:));
w = double(w(:));

%-------------------------------------------------

ms = 2*M + 1;
F = finufft1d1(t, w, +1, 1e-12, ms);
F = F(:);

a_vec = (1/N) * F;
inputssss = -M:M; 

%----------------------------------------

a = @(s) eval_series_finufft(s, a_vec);

toc;
tic;





%----------------------------------------------------------

Nint = 2^18;                % large power of two
sint = linspace(-pi,pi,Nint+1);
sint(end) = [];

a_vals = a(sint);

integrand_vals = 1./(2*pi*a_vals);

normalisation = (2*pi/Nint) * sum(integrand_vals);

a_normalised = @(s) normalisation .* a(s);
a = @(s) a_normalised(s);

toc;



%--------------------------------------------


function vals = eval_series_finufft(s, a_vec)

    xx = double(s(:));
    xx = mod(xx + pi, 2*pi) - pi;

    cc = complex(a_vec(:));

    opts = struct();  

    vals = finufft1d2(xx, +1, 1e-12, cc, opts);

    vals = reshape(vals, size(s));
end

end