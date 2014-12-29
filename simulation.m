num_observations = 500;
min_time = 0;
max_time = 10^9;
min_components = 1;
max_components = 6;
min_period = 0.1;
max_period = 50;
min_light = 12;
max_light = 18;
min_amp = 0;
max_amp = 2;
min_phi = 0;
max_phi = 2*pi;
err = 10^-5;
save = 0;
%for num_periods = 1:4
min_periods = 2;
max_periods = 2;

t = sort(unifrnd(min_time, max_time, num_observations, 1));
periodicities = floor(unifrnd(min_periods, max_periods+1, 1, 1));
p = sort(unifrnd(min_period, max_period, 1, periodicities), 'descend');
n = floor(unifrnd(min_components, max_components+1, 1, periodicities));
A = zeros(periodicities, max(n));
B = zeros(periodicities, max(n));
for p_i = 1:periodicities
    zs = zeros(1, max(n) - n(p_i));
    A(p_i, :) = [unifrnd(min_amp,max_amp,1,n(p_i))./10.^(1:n(p_i)) zs];
    B(p_i, :) = [unifrnd(min_amp,max_amp,1,n(p_i))./10.^(1:n(p_i)) zs];
end

A0 = unifrnd(min_light, max_light, 1, 1);
m = ones(num_observations, 1) * A0 + normrnd(0, err, num_observations, 1);
for ii = 1:num_observations
    for jj = 1:periodicities
        ph = t(ii) ./ p(jj);
        for kk = 1:n(jj)
            m(ii) = m(ii) + A(jj, kk) * sin(2*pi*kk*ph) ...
                          + B(jj, kk) * cos(2*pi*kk*ph);
        end
    end
end

figure(h)
set(gca,'YDir','reverse')
subplot(1, periodicities+1, 1)
plot(t, m, '+')
xlabel('Time')
ylabel('Magnitude')
title('Raw')
for jj = 1:periodicities
    subplot(1, periodicities+1, jj+1)
    x = mod(t ./ p(jj), 1);
    plot(x, m, '+')
    
    %[x, x_index] = sort(mod(t ./ p(d_i), 1));
    %ms = m(x_index);
    
    %[val, m_index] = max(m);
    %plot(mod(x+x(m_index),1), m, '+');
    
    %plot([x(m_index:length(m)); x(1:m_index-1)], ...
    %     [ms(m_index:length(m)); ms(1:m_index-1)], '+')
    
    hold on
    x = 0:.01:1;
    fx = A0;
    for kk = 1:n(jj)
        fx = fx + A(jj,kk) * sin(2*pi*kk*x) ...
                + B(jj,kk) * cos(2*pi*kk*x);
    end
    plot(x, fx, 'r', 'LineWidth', 2)
    hold off
    
    set(gca,'YTick',[])
    xlabel('Phase')
    title(['P = ' num2str(p(jj)) 'd'])
end
if save; saveas(h, ['plots/' int2str(num_periods) '.png']); end;
%end
