function [h,k,phi] = reconstruction_diff_re_large(a,number_of_modes_of_h_one_side,resolution,theta_prime)


numb_of_disc = 10*number_of_modes_of_h_one_side;

N = numb_of_disc; 

tic;h = @(s) -log(a(s));
      
L = 2*pi;       
grid = linspace(-pi, pi, N+1);
grid(end) = [];

samples = h(grid);

F_h = fft(samples)/N;


fourier_coefficients_of_h = F_h(1:number_of_modes_of_h_one_side+1);

inputsssss = 1:number_of_modes_of_h_one_side;

coeffs = fourier_coefficients_of_h(2:end);
n = inputsssss;

k = @(s) k_scalar_eval(s, fourier_coefficients_of_h(1), coeffs, n);

phi = 0;toc;



    function val = k_scalar_eval(s, c0, cp, n)

       
        if numel(s) ~= 1
            error('k expects scalar input.');
        end
    
        n  = n(:);      
        cp = cp(:);     
    
        z = exp(1i * s * n.');   
    
        val = c0 + 2 * (z * cp); 

    end
end


