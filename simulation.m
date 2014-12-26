num_observations = 500;
min_time = 0;
max_time = 10^6;
min_periods = 1;
max_periods = 3;
min_components = 1;
max_components = 7;
min_period = 0.1;
max_period = 50;
min_light = 12;
max_light = 18;
min_amp = -1;
max_amp = 1;
min_phi = 0;
max_phi = 2*pi;
err = 0;%10^-10;

t = sort(unifrnd(min_time, max_time, num_observations, 1));
periodicities = floor(unifrnd(min_periods, max_periods+1, 1, 1));
p = sort(unifrnd(min_period, max_period, 1, periodicities), 'descend');
n = floor(unifrnd(min_components, max_components+1, 1, periodicities));
A = zeros(periodicities, max(n), 2);
for p_i = 1:periodicities
    zs = zeros(1, max(n) - n(p_i));
    %A(p_i, :, 1) = [unifrnd(min_amp,max_amp,2,n(p_i)) ...
    %                    ./ 10 .^ [1:n(p_i); 1:n(p_i)] ...
    %                zeros(2, max(n) - n(p_i))]';%zs];
    A(p_i, :, 1) = [unifrnd(min_amp,max_amp,1,n(p_i))./10.^(1:n(p_i)) zs];
    A(p_i, :, 2) = [unifrnd(min_amp,max_amp,1,n(p_i))./10.^(1:n(p_i)) zs];
end

A0 = unifrnd(min_light, max_light, 1, 1);
m = ones(num_observations, 1) * A0 + normrnd(0, err, num_observations, 1);
for ii = 1:num_observations
    for jj = 1:periodicities
        ph = mod(t(ii) ./ p(jj), 1);
        for kk = 1:n(jj)
            for ll = 0:1
                m(ii) = m(ii) + A(jj,kk,ll+1) * sin(2*pi*kk*ph + ll*pi/2);
            end
            %m(ii) = m(ii) + ...
            %    A(jj, kk, 1) * cos(2*pi*kk * mod(t(ii) ./ p(jj), 1)) + ...
            %    A(jj, kk, 2) * sin(2*pi*kk * mod(t(ii) ./ p(jj), 1));
        end
    end
end

figure
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
        for ll = 0:1
            ph = x;%mod(x+x(m_index),1);
            fx = fx + A(jj,kk,ll+1) * sin(2*pi*kk*ph + ll*pi/2);
        end
    end
    plot(x, fx, 'r', 'LineWidth', 2)
    hold off
    
    set(gca,'YTick',[])
    xlabel('Phase')
    title(['P = ' num2str(p(jj)) 'd'])
end
