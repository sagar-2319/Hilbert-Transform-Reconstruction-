M = 200001;
thetaa = linspace(-pi,pi,M);
thetaa(end) = [];

N = M-1;

% ----------------------
pool = gcp('nocreate');
if isempty(pool)
    parpool(8);
end

numWorkers = gcp().NumWorkers;
numBlocks  = 4 * numWorkers;

edges = round(linspace(1, N+1, numBlocks+1));

% ---------------------
k_blocks = cell(numBlocks,1);

% ----------------------------
D = parallel.pool.DataQueue;
progressCount = 0;
startTime = tic;
lastPct = -1;

afterEach(D, @updateProgress);
fprintf('Progress:   0%% | ETA: --:--');

tic;

parfor b = 1:numBlocks

    s = edges(b);
    e = edges(b+1)-1;

    tmp = zeros(1, e-s+1);

    for j = 1:(e-s+1)
        tmp(j) = k_func(thetaa(s+j-1));
    end

    k_blocks{b} = tmp;

    send(D, e-s+1); 
end

%------------------------------
k_theta = zeros(1,N);

for b = 1:numBlocks
    s = edges(b);
    e = edges(b+1)-1;
    k_theta(s:e) = k_blocks{b};
end

toc;
fprintf('\n');

% ----------------------------

f_theta = exp(k_theta);

f_hat = fft(f_theta) / N;
neg_energy = norm(f_hat(floor(N/2)+2:end));

phi_hat = zeros(size(f_hat));
n = 0:N-1;
phi_hat(2:end) = f_hat(1:end-1) ./ (n(1:end-1) + 1);

Phi_theta = ifft(phi_hat * N)+.466353;

x1 = real(Phi_theta);
y1 = imag(Phi_theta);

x_plot = [x1, x1(1)];
y_plot = [y1, y1(1)];

% function data
gridd = linspace(0,2*pi, 10000);

samples = (a__.*cos(gridd) + 1i.*b__.*sin(gridd))+0.4;                       % function values

x2 = real(samples);
y2 = imag(samples);



% plot together
figure(30);
plot(x_plot, y_plot,'r', 'LineWidth', 1.5);%, x2, y2, '-');%,x3,y3,'.')
%plot(x2, y2, '-k','LineWidth   ', 1.5);
axis equal
grid on

ax = gca;
ax.XAxisLocation = 'origin';   % x-axis at y = 0
ax.YAxisLocation = 'origin';   % y-axis at x = 0

maxVal = max([abs(x2(:)); abs(y2(:))]);
maxVal = .8;
factor = 1;   % increase this (e.g., 1.5, 2) to zoom out more
xlim([-factor*maxVal factor*maxVal])
ylim([-factor*maxVal factor*maxVal])

ax.TickDir = 'in';             % ticks point inside
box off                        % optional (clean look)

legend('reconstructed', 'original')
xlabel('x')
ylabel('y')


%% ===== progress callback =====
function updateProgress(nDone)
    persistent count lastPct tStart Ntotal

    if isempty(count)
        count = 0;
        lastPct = -1;
        tStart = evalin('base','startTime');
        Ntotal = evalin('base','N');
    end

    count = count + nDone;
    pct = floor(100*count/Ntotal);

    if pct > lastPct
        elapsed = toc(tStart);
        rate = count / max(elapsed,eps);
        remaining = (Ntotal - count) / rate;

        % format ETA
        if remaining > 3600
            etaStr = sprintf('%dh %02dm', floor(remaining/3600), ...
                             floor(mod(remaining,3600)/60));
        elseif remaining > 60
            etaStr = sprintf('%dm %02ds', floor(remaining/60), ...
                             floor(mod(remaining,60)));
        else
            etaStr = sprintf('%ds', floor(remaining));
        end

        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
        fprintf('Progress: %3d%% | ETA: %s', pct, etaStr);

        lastPct = pct;
    end
end