%% Simulation
num_points = 20;
err = 0.1;
x_max = 10*2*pi;
t = 0:0.1:x_max;
p2 = 0.71;
m = sin(t) + 0.4 * sin(2*t+pi/4) ...
    + 0.2 * sin(p2 * t + pi) ...
    + 0.1 * sin(p2 * 8 * t + pi/15);

ust = sort(unifrnd(0, x_max, num_points, 1));
noise = unifrnd(0, err, num_points, 1);
usm = sin(ust) + 0.4*sin(2*ust+pi/4) ...
      + 0.2 * sin(p2 * ust + pi) ...
      + 0.1 * sin(p2 * 8 * ust+pi/15) ...
      + noise;

num_p1 = 4;
num_p2 = 16;
X = ones(length(ust), 1+num_p1+num_p2);
for kk = 2:2:num_p1
    X(:,kk)   = sin(kk/2*ust);
    X(:,kk+1) = cos(kk/2*ust);
end
for kk = 2:num_p2
    X(:,kk+num_p1)   = sin(p2 * kk/2*ust);
    X(:,kk+num_p1+1) = cos(p2 * kk/2*ust);
end

w = (X' * X) \ X' * usm;

wm = ones(1, length(t)) * w(1);
for kk = 2:2:num_p1
    wm = wm + w(kk)   * sin(kk/2*t) ...
            + w(kk+1) * cos(kk/2*t);
end
for kk = 2:2:num_p2
    wm = wm + w(kk+num_p1)   * sin(p2 * kk/2*t) ...
            + w(kk+num_p1+1) * cos(p2 * kk/2*t);
end

[w, FitInfo] = lasso(X, usm);%, 'cv', 10);
lm = ones(1, length(t)) * FitInfo.Intercept(1);
for kk = 2:2:num_p1
    lm = lm + w(kk)   * sin(kk/2*t) ...;
            + w(kk+1) * cos(kk/2*t);
end
for kk = 2:2:num_p2
    lm = lm + w(kk+num_p1)   * sin(p2*kk/2*t) ...
            + w(kk+num_p1+1) * cos(p2*kk/2*t);
end

ymax = max([max(m), max(wm), max(usm), max(lm)]) + 4*err;
ymin = min([min(m), min(wm), min(usm), min(lm)]) - 4*err;

% Plot time series
figure;
plot(t, m, 'color', [.5 0 0], 'LineWidth', 1.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
matlab2tikz('mpo.tikz', 'height', '\figureheight', ...
                         'width', '\figurewidth');

% Plot points along curve
figure;
plot(t, m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
hold on;
errorbar(ust, usm, ones(1,num_points)*2*err, '.', 'MarkerSize', 7.5,...
         'color', [.5 0 0]);
hold off;
matlab2tikz('mpo-points.tikz', 'height', '\figureheight', ...
                               'width',  '\figurewidth');

% Plot OLS
figure;
plot(t, m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1);
hold on;
plot(t, wm, '-', 'color', [0 0 0], 'LineWidth', 0.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
errorbar(ust, usm, ones(1,num_points)*2*err, '.', 'MarkerSize', 7.5,...
         'color', [.5 0 0]);
hold off;
matlab2tikz('mpo-badfit.tikz', 'height', '\figureheight', ...
                               'width',  '\figurewidth');

% Plot LASSO
figure;
plot(t, m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1);
hold on;
plot(t, lm, '-', 'color', [0 0 0], 'LineWidth', 0.5);
xlabel('t');
ylabel('m(t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca, 'box', 'off');
set(gca, 'ylim', [ymin, ymax]);
set(gca, 'xlim', [0, x_max]);
errorbar(ust, usm, ones(1,num_points)*2*err, '.', 'MarkerSize', 7.5,...
         'color', [.5 0 0]);
hold off;
matlab2tikz('mpo-lasso.tikz', 'height', '\figureheight', ...
                              'width', '\figurewidth');

% Plot in Fourier space
figure;
plot(sin(t), cos(t), '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1);
xlabel('\sin(\omega_1t)');
ylabel('\cos(\omega_1t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
hold on;
plot(X(:,2), X(:,3), '.', 'color', [.5 0 0]);
hold off;
matlab2tikz('mpo-fourier-space.tikz', 'height', '\figureheight', ...
                                      'width', '\figurewidth');

% Plot 3D view
%subplot(2, 2, 2);
figure;
plot3(sin(t), cos(t), m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 0.5);
xlabel('\sin(\omega_1t)');
ylabel('\cos(\omega_1t)');
zlabel('m(t)');
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'zticklabel',[]);
hold on;
scatter3(X(:,2), X(:,3), usm, 25, ...
         repmat([0.5,0,0],numel(usm),1), 'marker', '.');
for ii = 1:length(usm)
    plot3([sin(ust(ii)) sin(ust(ii))], [cos(ust(ii)) cos(ust(ii))], ...
        [usm(ii)-err usm(ii)+err], 'color', [0.5, 0, 0]);
end
grid on;
hold off;
matlab2tikz('mpo-fourier-scaled.tikz', 'height', '\figureheight', ...
                                       'width', '\figurewidth');
%{
% Pass plane through using LASSO coefficients 
figure;
plot3(sin(t), cos(t), m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 0.5);
xlabel('\sin(\omega_1t)');
ylabel('\cos(\omega_1t)');
zlabel('m(t)');
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'zticklabel',[]);
hold on;
scatter3(x(:,2), x(:,3), usm, 25, ...
         repmat([0.5,0,0],numel(usm),1), 'marker', '.');
for ii = 1:length(usm)
    plot3([sin(ust(ii)) sin(ust(ii))], [cos(ust(ii)) cos(ust(ii))], ...
        [usm(ii)-err usm(ii)+err], 'color', [0.5, 0, 0]);
end
%{
P1 = [-1, -1, FitInfo.Intercept(1) - w(2) - w(3)];
P2 = [ 1, -1, FitInfo.Intercept(1) + w(3) - w(3)];
P3 = [ 1,  1, FitInfo.Intercept(1) + w(3) + w(3)];
normal = cross(P1-P2, P1-P3);
syms asdf1 asdf2 asdf3
asdfP = [asdf1,asdf2,asdf3];
realdot = @(u, v) u*transpose(v);
planefunction = realdot(normal, asdfP-P1);
%}
%{
[sinx, cosx] = meshgrid(-10:0.1:10);
z = w(2)*sinx + w(3)*cosx + FitInfo.Intercept(1);
mesh(sinx, cosx, z);
set(gca, 'ylim', [-1 1]);
set(gca, 'xlim', [-1 1]);
set(gca, 'zlim', [-2 2]);
%}
grid on;
hold off;
%}
subplot(2, 2, 3);
plot(sin(t), sin(2*t), '--', 'color', [0.5 0.5 0.5], 'LineWidth', 1);
xlabel('\sin(\omega_1t)');
ylabel('\sin(2\omega_1t)');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
hold on;
plot(X(:,2), X(:,4), '.', 'color', [.5 0 0]);
hold off;

subplot(2, 2, 4);
plot3(sin(t), sin(2*t), m, '--', 'color', [0.5 0.5 0.5], 'LineWidth', 0.5);
xlabel('\sin(\omega_1t)');
ylabel('\sin(2\omega_1t)');
zlabel('m(t)');
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'zticklabel',[]);
hold on;
scatter3(X(:,2), X(:,4), usm, 15, ...
         repmat([0.5,0,0],numel(usm),1), 'marker', '.');
for ii = 1:length(usm)
    plot3([sin(ust(ii)) sin(ust(ii))], [sin(2*ust(ii)) sin(2*ust(ii))], ...
        [usm(ii)-err usm(ii)+err], 'color', [0.5, 0, 0]);
end
grid on;
hold off;
matlab2tikz('mpo-fourier.tikz', 'height', '\figureheight', ...
                                'width',  '\figurewidth');



                                  

%% OGLE LMC CEP 0227
data = fscanf(fopen('OGLE-LMC-CEP-0227.dat','r'), '%f %f %f', [3 Inf])';
t = data(:, 1);
m = data(:, 2);
e = data(:, 3);
p = [309.404 3.797086];
%p = [309.669 3.797086];
%{
figure
subplot(1, length(p)+1, 1)
errorbar(t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [min(t) max(t)])
set(gca, 'ylim', [min(m)-2*max(e) max(m)+2*max(e)])
xlabel('t')
ylabel('m(t)')
title('OGLE-LMC-CEP-0227')
for jj = 1:length(p)
    subplot(1, length(p)+1, jj+1)
    ph = mod(t ./ p(jj), 1);
    errorbar(ph, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
    set(gca, 'YDir', 'reverse')
    set(gca, 'xlim', [0 1])
    set(gca, 'ylim', [min(m)-2*max(e) max(m)+2*max(e)])
    set(gca, 'YTick', [])
    xlabel('Phase')
    title(['P = ' num2str(p(jj)) 'd'])
end
%}
num_p1 = 100;
num_p2 = 12;
num_fit = 40;
X = ones(length(m), 1+num_p1+num_p2);
for kk = 2:2:num_p1
    X(:,kk)   = sin(2*pi*kk/2*mod(t/p(1),1));
    X(:,kk+1) = cos(2*pi*kk/2*mod(t/p(1),1));
end
for kk = 2:num_p2
    X(:,kk+num_p1)   = sin(2*pi*kk/2*mod(t/p(2),1));
    X(:,kk+num_p1+1) = cos(2*pi*kk/2*mod(t/p(2),1));
end

[w, FitInfo] = lasso(X, m, 'CV', 3);
w = w(:,num_fit);
%w = (x' * x) \ x' * m;

figure
subplot(1, length(p)+1, 1)
ph = min(t):1:max(t);
lm = ones(1, length(ph)) * FitInfo.Intercept(num_fit);
for kk = 2:2:num_p1
    lm = lm + w(kk)   * sin(2*pi*kk/2*mod(ph/p(1),1)) ...;
            + w(kk+1) * cos(2*pi*kk/2*mod(ph/p(1),1));
end
for kk = 2:2:num_p2
    lm = lm + w(kk+num_p1)   * sin(2*pi*kk/2*mod(ph/p(2),1)) ...
            + w(kk+num_p1+1) * cos(2*pi*kk/2*mod(ph/p(2),1));
end
plot(ph, lm, 'LineWidth', 1, 'Color', [0 0 0])
hold on
errorbar(t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [min(t) max(t)])
set(gca, 'ylim', [min(m)-2*max(e) max(m)+2*max(e)])
xlabel('t')
ylabel('m(t)')
title('OGLE-LMC-CEP-0227 Raw Photometry')

subplot(1, length(p)+1, 2)
ph = mod(t ./ p(1), 1);
errorbar(ph, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [0 1])
set(gca, 'ylim', [min(m)-2*max(e) max(m)+2*max(e)])
hold on
ph = 0:.01:1;
lm = ones(1, length(ph)) * FitInfo.Intercept(num_fit);
for kk = 2:2:num_p1
    lm = lm + w(kk)   * sin(2*pi*kk/2*ph) ...
            + w(kk+1) * cos(2*pi*kk/2*ph);
end
plot(ph, lm, 'LineWidth', 0.5, 'Color', [0 0 0])
hold off
set(gca, 'YTick', [])
xlabel('Phase')
title(['P = ' num2str(p(1)) 'd'])

subplot(1, length(p)+1, 3)
ph = mod(t ./ p(2), 1);
errorbar(ph, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [0 1])
set(gca, 'ylim', [min(m)-2*max(e) max(m)+2*max(e)])
hold on
ph = 0:.01:1;
lm = ones(1, length(ph)) * FitInfo.Intercept(num_fit);
for kk = 2:2:num_p2
    lm = lm + w(kk+num_p1)   * sin(2*pi*kk/2*ph) ...
            + w(kk+num_p1+1) * cos(2*pi*kk/2*ph);
end
plot(ph, lm, 'LineWidth', 0.5, 'Color', [0 0 0])
hold off
set(gca, 'YTick', [])
xlabel('Phase')
title(['P = ' num2str(p(2)) 'd'])









%% PL
asdf = [''];
for ii = 1:101
    asdf = ['%f ' asdf];
end
ogle_i = fscanf(fopen('ceps.dat','r'), asdf, [101 Inf])';

smallp = log10(ogle_i(:,1)) < 1;
largep = ~smallp;
r2max = 0;
maxin = 0;
r2min = 1;
minin = 0;
for ii = 2:101
    mdl = fitlm(log10(ogle_i(:,1)), ogle_i(:,ii));
    mdls = fitlm(log10(ogle_i(smallp,1)), ogle_i(smallp,ii));
    mdll = fitlm(log10(ogle_i(~smallp,1)), ogle_i(~smallp,ii));
    [ii mdl.Rsquared.Ordinary ...
        mdls.Rsquared.Ordinary ...
        mdll.Rsquared.Ordinary]
    %{
    if mdl.Rsquared.Ordinary > r2max
        r2max = mdl.Rsquared.Ordinary;
        maxin = ii;
    end
    
    if mdl.Rsquared.Ordinary < r2min
        r2min = mdl.Rsquared.Ordinary;
        minin = ii;
    end
    %}
end

figure
scatter(log10(ogle_i(:,1)), ogle_i(:,2), 23, [.5 0 0], '.')
scatter(log10(ogle_i(:,1)), ogle_i(:,85), 23, [.5 0 0], '.')
set(gca, 'YDir', 'reverse')
xlabel('Log Period')
ylabel('m(\phi = 0)')







%% OGLE-SMC-LPV-12337
name = 'OGLE-SMC-LPV-12337';
p = [262.9];%[12.9751466];
data = fscanf(fopen([name '.dat'],'r'), '%f %f %f', [3 Inf])';
t = data(:, 1);
m = data(:, 2);
e = data(:, 3);

num_p1 = 16;
num_fit = 10;
X = ones(length(m), 1+num_p1);
for kk = 2:2:num_p1
    X(:,kk)   = sin(2*pi*kk/2*mod(t/p(1),1));
    X(:,kk+1) = cos(2*pi*kk/2*mod(t/p(1),1));
end
%[w, FitInfo] = lasso(x, m);%, 'CV', 10);
%w = w(:,num_fit);
w = (X' * X) \ X' * m;

figure
%subplot(1, length(p)+1, 2)
ph = mod(t ./ p(1), 1);
errorbar(ph, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [0 1])
set(gca, 'ylim', [min(m)-2*max(e) max(m)+2*max(e)])
hold on
ph = 0:.01:1;
lm = ones(1, length(ph)) * FitInfo.Intercept(num_fit);
for kk = 2:2:num_p1
    lm = lm + w(kk)   * sin(2*pi*kk/2*ph) ...
            + w(kk+1) * cos(2*pi*kk/2*ph);
end
plot(ph, lm, 'LineWidth', 0.5, 'Color', [0 0 0])
hold off
set(gca, 'ylim', [min(m)-max(e) max(m)+max(e)])
%set(gca, 'YTick', [])
ylabel('m(\phi)')
xlabel(['\phi (P = ' num2str(p(1)) 'd)'])
%title()
matlab2tikz([name '.tikz'], 'height', '\figureheight', ...
                            'width', '\figurewidth');

figure
%subplot(1, length(p)+1, 1)
%{
ph = min(t):1:max(t);
lm = ones(1, length(ph)) * FitInfo.Intercept(num_fit);
for kk = 2:2:num_p1
    lm = lm + w(kk)   * sin(2*pi*kk/2*mod(ph/p(1),1)) ...;
            + w(kk+1) * cos(2*pi*kk/2*mod(ph/p(1),1));
end
plot(ph, lm, 'LineWidth', 0.1, 'Color', [0 0 0])
hold on
%}
errorbar(t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1);
%hold off
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [min(t)-p(1) max(t)+p(1)])
set(gca, 'ylim', [min(m)-max(e) max(m)+max(e)])
xlabel('t')
ylabel('m(t)')
%title('OGLE-BLG-CEP-31')
matlab2tikz([name '-photometry.tikz'], 'height', '\figureheight', ...
                                      'width', '\figurewidth');










%% OGLE-LMC-CEP-0209.dat
% Inputs
name = 'OGLE-LMC-CEP-0209';
p = 3.1227238;
num_p1 = 32;
num_repeats = 2;
resolution = .001;-

% Parse data
data = fscanf(fopen([name '.dat'],'r'), '%f %f %f', [3 Inf])';
t = data(:, 1);
m = data(:, 2);
e = data(:, 3);

% Build feature matrix
X = ones(length(m), 1+num_p1);
for kk = 2:2:num_p1
    X(:,kk)   = sin(2*pi*kk/2*mod(t/p(1),1));
    X(:,kk+1) = cos(2*pi*kk/2*mod(t/p(1),1));
end

% Fit OLS and Lasso
ph = 0:resolution:num_repeats;
[w, FitInfo] = lasso(X, m, 'CV', 3);
w = w(:,FitInfo.IndexMinMSE);
lm = ones(1, length(ph)) * FitInfo.Intercept(FitInfo.IndexMinMSE);
for kk = 2:2:num_p1
    lm = lm + w(kk)   * sin(2*pi*kk/2*ph) ...
            + w(kk+1) * cos(2*pi*kk/2*ph);
end

w = (X' * X) \ X' * m;
ols = ones(1, length(ph)) * w(1);
for kk = 2:2:num_p1
    ols = ols + w(kk)   * sin(2*pi*kk/2*ph) ...
              + w(kk+1) * cos(2*pi*kk/2*ph);
end

ymax = max([max(m), max(lm), max(ols)]) + 4*max(e);
ymin = min([min(m), min(lm), min(ols)]) - 4*max(e);

% Plot unphased data
figure
errorbar(t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1);
set(gca, 'YDir', 'reverse')
offset = .05 * (max(t)-min(t));
set(gca, 'xlim', [min(t)-offset max(t)+offset])
set(gca, 'ylim', [ymin ymax])
xlabel('t')
ylabel('m(t)')
matlab2tikz([name '-photometry.tikz'], 'height', '\figureheight', ...
                                      'width', '\figurewidth');

% Plot phased data
figure
ph = mod(t ./ p(1), num_repeats);
errorbar(ph, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [0 num_repeats])
set(gca, 'ylim', [ymin ymax])
ylabel('m(\phi)')
xlabel(['\phi (P = ' num2str(p(1)) 'd)'])
matlab2tikz([name '-phased.tikz'], 'height', '\figureheight', ...
                                   'width', '\figurewidth');

% Plot OLS
figure
ph = mod(t ./ p(1), num_repeats);
errorbar(ph, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [0 num_repeats])
hold on
ph = 0:resolution:num_repeats;
plot(ph, ols, 'LineWidth', 0.5, 'Color', [0 0 0])
hold off
set(gca, 'ylim', [ymin ymax])
ylabel('m(\phi)')
xlabel(['\phi (P = ' num2str(p(1)) 'd)'])
matlab2tikz([name '-ols.tikz'], 'height', '\figureheight', ...
                            'width', '\figurewidth');

% Plot lasso
figure
ph = mod(t ./ p(1), num_repeats);
errorbar(ph, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1)
set(gca, 'YDir', 'reverse')
set(gca, 'xlim', [0 num_repeats])
hold on
ph = 0:resolution:num_repeats;
plot(ph, lm, 'LineWidth', 0.5, 'Color', [0 0 0])
hold off
set(gca, 'ylim', [ymin ymax])
ylabel('m(\phi)')
xlabel(['\phi (P = ' num2str(p(1)) 'd)'])
matlab2tikz([name '-lasso.tikz'], 'height', '\figureheight', ...
                            'width', '\figurewidth');







%% OGLE-LMC-CEP-0207.dat
% Inputs
name = 'OGLE-LMC-CEP-0432';
p = [3.8842199 2.7931983];
max_A = 10;
num_repeats = 2;
resolution = 1;

% Parse data
data = fscanf(fopen([name '.dat'], 'r'), '%f %f %f', [3 Inf])';
t = data(:, 1);
m = data(:, 2);
e = data(:, 3);
ph = [(min(t):resolution:max(t))'; t];
ph_s = (0:.1:ceil(prod(p)*4))';

% Build feature matrix
ks = cell(1,numel(p)); 
[ks{:}] = ndgrid(0:max_A);
ks = reshape(cat(numel(p)+1, ks{:}), [], numel(p));
X   = [cos(t*(ks*(2*pi./p'))') sin(t*(ks(2:length(ks),:)*(2*pi./p'))')];
Xph = [cos(ph*(ks*(2*pi./p'))') sin(ph*(ks(2:length(ks),:)*(2*pi./p'))')];
Xs  = [cos(ph_s*(ks*(2*pi./p'))') sin(ph_s*(ks(2:length(ks),:)*(2*pi./p'))')];

% Fit OLS and Lasso
[lw, FitInfo] = lasso(X, m, 'CV', 3);
lm = Xph*lw(:,FitInfo.IndexMinMSE) + FitInfo.Intercept(FitInfo.IndexMinMSE);
lms = Xs*lw(:,FitInfo.IndexMinMSE) + FitInfo.Intercept(FitInfo.IndexMinMSE);
%lsw = (X' * X) \ X' * m;
%lsm = Xph*lsw;
ymax = max([max(m), max(lm)]) + 4*max(e);%, max(lsm)]) + 4*max(e);
ymin = min([min(m), min(lm)]) - 4*max(e);%, min(lsm)]) - 4*max(e);

% Plot unphased data
figure
subplot(1, 2, 1)
errorbar(t, m, e, '.', 'color', [.5 0 0], 'MarkerSize', 1);
set(gca, 'YDir', 'reverse')
offset = .05 * (max(t)-min(t));
set(gca, 'xlim', [min(t)-offset max(t)+offset])
set(gca, 'ylim', [ymin ymax])
xlabel('t')
ylabel('m(t)')
subplot(1, 2, 2)
plot(ph_s, lms, 'LineWidth', 0.5, 'Color', [0 0 0])
set(gca, 'xlim', [0, max(ph_s)]);
set(gca, 'ylim', [ymin ymax]);
ylabel('m(\phi)')
xlabel(['\phi'])
matlab2tikz([name '-photometry.tikz'], 'height', '\figureheight', ...
                                      'width', '\figurewidth');

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
matlab2tikz([name '-ols.tikz'], 'height', '\figureheight', ...
                            'width', '\figurewidth');
