function [r2 lw] = lightcurve( name, p, max_A, params )

if isfield(params, 'num_repeats')
    num_repeats = params.num_repeats;
else
    num_repeats = 2;
end

if isfield(params, 'phase_resolution')
    phase_resolution = params.phase_resolution;
else
    phase_resolution = 1;
end

if isfield(params, 'ts_resolution')
    ts_resolution = params.ts_resolution;
else
    ts_resolution = 1;
end

if isfield(params, 'windows')
    windows = params.windows;
else
    windows = 1;
end

if isfield(params, 'save_plot')
    save_plot = params.save_plot;
else
    save_plot = false;
end

% Parse data
data = fscanf(fopen([name '.dat'], 'r'), '%f %f %f', [3 Inf])';
t = data(:, 1);
m = data(:, 2);
e = data(:, 3);
ph = [(min(t) : phase_resolution : max(t))'; t];
ts = (floor(min(t)) : ts_resolution : ...
      ceil(min(t)+prod(windows*p)*windows))';

% Find wavenumbers (k) and angular frequencies (omega)
ks = cell(1,numel(p)); 
[ks{:}] = ndgrid(0:max_A);
ks = reshape(cat(numel(p)+1, ks{:}), [], numel(p));
omega = (2*pi./p');

% Build feature matrix (X)
cos_freqs = (ks*omega)';
sin_freqs = (ks(2:length(ks),:)*omega)'; % remove column of 0s
X   = [cos(t  * cos_freqs) sin(t  * sin_freqs)];
Xph = [cos(ph * cos_freqs) sin(ph * sin_freqs)];
Xts = [cos(ts * cos_freqs) sin(ts * sin_freqs)];

% Fit OLS and Lasso
[lw, FitInfo] = lasso(X, m, 'CV', 3);
lw = [FitInfo.Intercept(FitInfo.IndexMinMSE); ...
       lw(2:length(lw), FitInfo.IndexMinMSE)];
lm = Xph*lw;
%lsw = (X' * X) \ X' * m;
%lsm = Xph*lsw;
ymax = max([max(m), max(lm)]) + 4*max(e);%, max(lsm)]) + 4*max(e);
ymin = min([min(m), min(lm)]) - 4*max(e);%, min(lsm)]) - 4*max(e);

SStot = var(m) * length(m);
SSres = sum((X*lw - m).^2);
r2 = 1 - SSres/SStot;

% Plot unphased data
figure
subplot(1, 2, 1)
errorbar(t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1);
offset = .05 * (max(t)-min(t));
set(gca, 'xlim', [min(t)-offset max(t)+offset])
set(gca, 'ylim', [ymin ymax])
set(gca, 'YDir', 'reverse')
xlabel('t')
ylabel('m(t)')
subplot(1, 2, 2)
plot(ts, Xts*lw, 'LineWidth', 0.5, 'Color', [0 0 0])
set(gca, 'xlim', [min(ts) max(ts)]);
set(gca, 'ylim', [ymin ymax]);
set(gca, 'YDir', 'reverse')
xlabel('t')
ylabel('m(t)')
if save_plot
    matlab2tikz([name '-photometry.tikz'], 'height', '\figureheight', ...
                                           'width',  '\figurewidth');
end

% Fit!
figure
for p_i = 1:length(p)
    subplot(1, length(p), p_i)
    phased = sortrows([mod(ph ./ p(p_i), num_repeats) lm], 1);
    plot(phased(:,1), phased(:,2), ...
         'LineWidth', 0.5, 'Color', [0 0 0])
    set(gca, 'YDir', 'reverse')
    set(gca, 'xlim', [0 num_repeats])
    set(gca, 'ylim', [ymin ymax])
    ylabel('m(\phi)')
    xlabel(['\phi (P = ' num2str(p(p_i)) 'd)'])
    hold on
    ph_t = mod(t ./ p(p_i), num_repeats);
    errorbar(ph_t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
    hold off
end
if save_plot
    matlab2tikz([name '-phased.tikz'], 'height', '\figureheight', ...
                                       'width',  '\figurewidth');
end

end
