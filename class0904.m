function class0904
A = 1;
C = 1;

N=100;
sigma_x0 = 150;
x0 = mvnrnd(0, sigma_x0, 1);
x_minus = [x0];
E_minus = [1];

% сгенерируем помехи
sigma_v = 30;
sigma_w = 20;

v = zeros(1, N);
for i = 1 : N
    v(i) = mvnrnd(0, sigma_v);
end

w = zeros(1, N);
for i = 1 : N
    w(i) = mvnrnd(0, sigma_w);
end

% x и y
x  = zeros(1, N);
y  = zeros(1, N);
x(1) = x0;

y(1) = C * x(1) + w(1);
for i = 2 : N
    x(i) = A * x(i - 1) + v(i - 1);
    y(i) = C * x(i) + w(i);
end

% строим фильтр Каллмана
for k=2:N
    x_minus = [x_minus; x_minus(k-1) + (E_minus(k-1)*(1/(sigma_w+E_minus(k-1)))) * (y(k-1)-x_minus(k-1))];
    E_minus = [E_minus; 1/(1/E_minus(k-1) + 1/20) + 30];
end

% стационарный фильтр Каллмана
E_star = 15 + 5 * sqrt(33);
x_minus_star = [x0];
K_k = E_star / (E_star + sigma_w);
for k=2:N
    x_minus_star = [x_minus_star; x_minus_star(k-1) + K_k * (y(k-1)-x_minus_star(k-1))];
end

% x и y
x_bezpomeh  = zeros(1, N);
y_bezpomeh  = zeros(1, N);
x_bezpomeh(1) = x0;

y_bezpomeh(1) = C * x_bezpomeh(1);
for i = 2 : N
    x_bezpomeh(i) = A * x_bezpomeh(i - 1);
    y_bezpomeh(i) = C * x_bezpomeh(i);
end

% фильтр Каллмана
x_minus_bezpomeh = [x0];
E_minus_bezpomeh = [1];
for k=2:N
    x_minus_bezpomeh = [x_minus_bezpomeh; x_minus_bezpomeh(k-1) + y_bezpomeh(k-1)-x_minus_bezpomeh(k-1)];
end
disp(x0)

% апостериорная оценка
x_plus = [x0];
E_plus =[sigma_w/(sigma_w+1)];

for k=2:N
    x_plus = [x_plus; x_plus(k-1) + (E_plus(k-1) + sigma_v)/(sigma_w + E_plus(k-1) + sigma_v) * (y(k)-x_plus(k-1))];
    E_plus = [E_plus; inv( 1/(E_plus(k-1) + sigma_v) + 1/20)];
end
figure;
subplot(1, 1, 1)
    plot(1:N, x, 'g', 1:N, x_minus, 'b', 1:N, x_minus_star, 'r')
    grid on;
    xlabel('k');
    ylabel('x');
    legend({'Истинное x', 'Априорная оценка x', 'Априорн. оценка через стационарный фильтр'}, 'Location', 'northeast');
figure;
    subplot(1, 1, 1)
    plot(1:N, x_minus_star-x_minus, 'b');
    legend('Невязка для априорной оценки и апр оценки через стационарный фильтр', 'Location', 'northeast');
figure;
    subplot(1, 1, 1)
    plot(1:N, x_bezpomeh, 'b', 1:N, x_minus_bezpomeh, 'r');
    legend({'Истинное x без помех', 'Априорная оценка x без помех'}, 'Location', 'northeast');
figure;
subplot(1, 1, 1)
    plot(1:N, x, 'g', 1:N, x_minus, 'b', 1:N, x_plus, 'black')
    grid on;
    xlabel('k');
    ylabel('x');
    legend({'Истинное x', 'Априорная оценка x', 'Апостериорная оценка'}, 'Location', 'northeast');
end