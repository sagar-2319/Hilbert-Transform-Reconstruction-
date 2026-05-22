function [L_0_integral,L_p_integral,L_n_integral,E_0_func,E_p_func,E_n_func] = E_diff_re(number_of_modes_of_L_one_side,L_0,L_p,L_n)
%E_DIFF_R Summary of this function goes here

inputsss = -number_of_modes_of_L_one_side:number_of_modes_of_L_one_side;

inverse_n = inputsss;
inverse_n(idx(0,number_of_modes_of_L_one_side)) = 1;
inverse_n = 1./inverse_n;

L_0_integral_coefficient = L_0;
L_0_integral_coefficient(idx(0,number_of_modes_of_L_one_side)) = 0;
L_0_integral_coefficient = L_0_integral_coefficient.*inverse_n;

L_p_integral_coefficient = L_p;
L_p_integral_coefficient(idx(0,number_of_modes_of_L_one_side)) = 0;
L_p_integral_coefficient = L_p_integral_coefficient.*inverse_n;

L_n_integral_coefficient = L_n;
L_n_integral_coefficient(idx(0,number_of_modes_of_L_one_side)) = 0;
L_n_integral_coefficient = L_n_integral_coefficient.*inverse_n;



L_0_integral = @(s) arrayfun(@(t) ( (L_0(idx(0,number_of_modes_of_L_one_side))).*t -1i.*((exp(1i.*inputsss.*t)-1) * L_0_integral_coefficient.') ...
    ) .* (t >= -pi/3 & t < pi/3), s);
L_p_integral = @(s) arrayfun(@(t) ( (L_p(idx(0,number_of_modes_of_L_one_side))).*(t-2*pi/3) -1i.*((exp(1i.*inputsss.*t) - exp(2.*pi.*1i.*inputsss./3)) ...
            * L_p_integral_coefficient.')) .* (t >= pi/3 & t <= pi), s);
L_n_integral = @(s) arrayfun(@(t) ...
    ( (L_n(idx(0,number_of_modes_of_L_one_side))).*(t+2*pi/3) -1i.*((exp(1i.*inputsss.*t) ...
             - exp(-2.*pi.*1i.*inputsss./3))  * L_n_integral_coefficient.') ) .* (t >= -pi & t < -pi/3), s);

E_0_func = @(s) exp(1/2.*L_0_integral(s)).*(s >= -pi/3 & s < pi/3);
E_p_func = @(s) exp(-1/2.*L_p_integral(s)).*(s >= pi/3 & s <= pi);
E_n_func = @(s) exp(-1/2.*L_n_integral(s)).*(s >= -pi & s < -pi/3);

