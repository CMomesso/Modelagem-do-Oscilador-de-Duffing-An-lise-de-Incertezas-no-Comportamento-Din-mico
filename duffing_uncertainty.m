clc; clear; close all;

% Tempo de simulação
t0 = 0.0;
t1 = 20.0;
Ndt = 400;
tspan = linspace(t0, t1, Ndt);

% Parâmetros fixos
delta = 0.3;     % amortecimento
alpha = -1.0;    % rigidez linear
gamma = 0;       % sem força externa
x0 = 0.5;        % condição inicial de posição
v0 = 0.0;        % condição inicial de velocidade
IC = [x0; v0];

% Monte Carlo
rng(2024);                   % semente aleatória para reprodutibilidade
Ns = 200;                    % número de amostras
lambda = 1.0;                % taxa da distribuição exponencial
beta_samples = exprnd(1/lambda, [Ns, 1]);  % distribuição exponencial

% Soluções
X = zeros(Ndt, Ns);

for i = 1:Ns
    beta = beta_samples(i);
    duffing = @(t, y) [y(2); -delta*y(2) - alpha*y(1) - beta*y(1)^3];
    [~, Y] = ode45(duffing, tspan, IC);
    X(:, i) = Y(:, 1);
end

% Estatísticas
X_mean = mean(X, 2);
X_std = std(X, 0, 2);
X_p95 = prctile(X, 97.5, 2);
X_p05 = prctile(X, 2.5, 2);

% Gráfico
figure;
fill([tspan fliplr(tspan)], [X_p95' fliplr(X_p05')], [0.9 0.9 0.9], 'EdgeColor', 'none'); hold on;
plot(tspan, X_mean, 'b', 'LineWidth', 2);
plot(tspan, X_mean + X_std, 'r--', 'LineWidth', 1.5);
plot(tspan, X_mean - X_std, 'r--', 'LineWidth', 1.5);
xlabel('Tempo (s)');
ylabel('Deslocamento x(t)');
title('Oscilador de Duffing com Incerteza Exponencial em \beta');
legend('95% intervalo', 'Média', '±1 desvio padrão');
grid on;
