function [r2 lw] = lightcurve( name, p, max_A, params )

if isfield(params, 'num_repeats')
    num_repeats = params.num_repeats;
else
    num_repeats = 2;
end

if isfield(params, 'phase_resolution')
    phase_resolution = params.phase_resolution;
else
    phase_resolution = 2;
end

if isfield(params, 'ts_resolution')
    ts_resolution = params.ts_resolution;
else
    ts_resolution = 0.1;
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
%ph = [(min(t) : phase_resolution : max(t))'; t];

% Find a nice place to plot the entire light curve
epochs = find(diff(t(t>2300)) > median(diff(t)) + mad(diff(t)));
epochs = epochs(epochs > length(t)*.05);
min_ts = floor(min(t));
max_ts = ceil(t(min(epochs(1),50)));
points_i = t > min_ts & t < max_ts;
ts = sort([(min_ts : ts_resolution : max_ts)'; t(points_i)]);

% Find wavenumbers (k) and angular frequencies (omega)
ks = cell(1,numel(p)); 
[ks{:}] = ndgrid(0:max_A);
ks = reshape(cat(numel(p)+1, ks{:}), [], numel(p));
omega = (2*pi./p');

% Build feature matrix (X)
cos_freqs = (ks*omega)';
sin_freqs = (ks(2:length(ks),:)*omega)'; % remove column of 0s
X   = [cos(t  * cos_freqs) sin(t  * sin_freqs)];
Xts = [cos(ts * cos_freqs) sin(ts * sin_freqs)];
%Xph = [cos(ph * cos_freqs) sin(ph * sin_freqs)];

% Fit OLS and Lasso
[lw, FitInfo] = lasso(X, m, 'CV', 3);
lw = [FitInfo.Intercept(FitInfo.IndexMinMSE); ...
       lw(2:size(lw, 1), FitInfo.IndexMinMSE)];
lm = Xts*lw;
%lm = Xph*lw;
%lsw = (X' * X) \ X' * m;
%lsm = Xph*lsw;
ymax = max([max(m), max(lm)]) + 1.5*max(e);%, max(lsm)]) + 4*max(e);
ymin = min([min(m), min(lm)]) - 1.5*max(e);%, min(lsm)]) - 4*max(e);

SStot = var(m) * length(m);
SSres = sum((X*lw - m).^2);
r2 = 1 - SSres/SStot;
%{
% Plot unphased data
figure
h(1) = subplot(2, 1, 1);
errorbar(t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1);
offset = .02 * (max(t)-min(t));
set(gca, 'xlim', [min(t)-offset max(t)+offset])
set(gca, 'ylim', [ymin ymax])
set(gca, 'YDir', 'reverse')
ylabel('m(t)')
h(2) = subplot(2, 1, 2);
plot(ts, Xts*lw, 'LineWidth', 0.5, 'Color', [0 0 0])
hold on
errorbar(t(points_i), m(points_i), e(points_i), ...
         '.', 'color', [.5 0 0], 'MarkerSize', 1);
hold off
set(gca, 'xlim', [min(ts) max(ts)])
xlabel('t (HJD-2450000)')
set(gca, 'YDir', 'reverse')
%set(gca, 'ytick', []);
ylabel('$\hat{\text{m}}$(t)')
linkaxes(h, 'y');
if save_plot
    matlab2tikz([name '-photometry.tikz'], 'height', '\figureheight', ...
                                           'width',  '0.8\figurewidth');
end

% Fit!
figure
for p_i = 1:length(p)
    subplot(length(p), 1, p_i)
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
%}
% Phase plots and time series
figure
for p_i = 1:length(p)
    h(p_i) = subplot(2, length(p), p_i);
    ph_t = mod(t ./ p(p_i), num_repeats);
    errorbar(ph_t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
    set(gca, 'YDir', 'reverse')
    set(gca, 'xlim', [0 num_repeats])
    set(gca, 'ylim', [ymin ymax])
    if p_i == 1
        ylabel('m(\phi)')
    else
        set(gca, 'yticklabel', [])
    end
    xlabel(['\phi (P = ' num2str(p(p_i)) 'd)'])
    set(gca,'xaxisLocation','top')
    set(gca,'Xtick',0:0.5:2)
    %set(gca,'XtickLabel',time(1:3:end))
end
h(length(p)+1) = subplot(2, length(p), length(p)+1:length(p)*2);
plot(ts, Xts*lw, 'LineWidth', 0.5, 'Color', [0 0 0])
hold on
errorbar(t(points_i), m(points_i), e(points_i), ...
         '.', 'color', [.5 0 0], 'MarkerSize', 1);
hold off
set(gca, 'xlim', [min(ts) max(ts)])
xlabel('t (HJD-2450000)')
set(gca, 'YDir', 'reverse')
%set(gca, 'ytick', []);
ylabel('$\hat{\text{m}}$(t)')
linkaxes(h, 'y');
if save_plot
    matlab2tikz([name '-phased-ts.tikz'], 'height', '\figureheight', ...
                                          'width',  '\figurewidth');
end

end
