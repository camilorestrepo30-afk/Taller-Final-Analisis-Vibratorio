%% NF1 - Barrido del coeficiente de amortiguamiento y tiempo de asentamiento (1%)
% m*x'' + c*x' + k*x = 0 ,  x(0)=1 m , v(0)=0
clear; clc; close all;

%% ---------- GRAFICAS CON FONDO BLANCO ----------
set(groot, 'defaultFigureColor', 'w');          % fondo de la ventana
set(groot, 'defaultAxesColor', 'w');            % fondo de los ejes
set(groot, 'defaultAxesXColor', 'k', 'defaultAxesYColor', 'k', 'defaultAxesZColor', 'k');
set(groot, 'defaultTextColor', 'k', 'defaultAxesGridColor', [0.15 0.15 0.15]);
try, set(groot, 'defaultLegendColor', 'w', 'defaultLegendTextColor', 'k'); end
set(groot, 'defaultAxesFontSize', 11);
try, set(groot, 'defaultFigureTheme', 'light'); end   % MATLAB R2025a+ (modo oscuro)

%% ---------- PARAMETROS (modificables) ----------
m   = 2.5;          % masa [kg]
k   = 100;          % rigidez [N/m]
cv  = 0:0.5:60;     % vector de coeficientes de amortiguamiento [N s/m]
x0  = 1;  v0 = 0;   % condiciones iniciales
tol = 0.01;         % criterio: 1% de x0
tf  = 100;          % tiempo maximo de simulacion [s]
dt  = 1e-3;         % paso de salida [s]

%% ---------- VALORES TEORICOS ----------
wn = sqrt(k/m);            % frecuencia natural [rad/s]
cc = 2*sqrt(k*m);          % amortiguamiento critico [N s/m]
fprintf('wn = %.4f rad/s ,  cc = %.4f N s/m\n', wn, cc);

%% ---------- BARRIDO ----------
t    = 0:dt:tf;
ts   = nan(size(cv));      % tiempo de asentamiento para cada c
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

for i = 1:length(cv)
    c = cv(i);
    f = @(t,y) [y(2); -(c*y(2) + k*y(1))/m];     % sistema de 1er orden
    [~,y] = ode45(f, t, [x0; v0], opts);
    x = y(:,1);

    fuera = find(abs(x) >= tol*x0);               % muestras fuera de la banda
    if isempty(fuera)
        ts(i) = 0;
    elseif fuera(end) == length(t)
        ts(i) = NaN;                              % nunca se asienta (c = 0)
    else
        j  = fuera(end);                          % ultima muestra fuera de la banda
        % interpolacion lineal para el instante exacto en que |x| = 1%
        a1 = abs(x(j)); a2 = abs(x(j+1));
        ts(i) = t(j) + (a1 - tol*x0)/(a1 - a2)*dt;
    end
end

%% ---------- RESULTADOS ----------
[tsmin, imin] = min(ts);
fprintf('Tiempo minimo de asentamiento: %.4f s con c = %.1f N s/m (zeta = %.3f)\n', ...
        tsmin, cv(imin), cv(imin)/cc);
ts_cc = interp1(cv, ts, cc);
fprintf('En el punto critico c = cc = %.3f N s/m  ->  ts = %.4f s\n', cc, ts_cc);

figure;
plot(cv, ts, 'b.-'); hold on; grid on;
plot([cc cc], [0 max(ts)], 'r--', 'LineWidth', 1.5);
plot(cc, ts_cc, 'ro', 'MarkerFaceColor', 'r');
plot(cv(imin), tsmin, 'ks', 'MarkerFaceColor', 'g');
xlabel('c [N s/m]'); ylabel('t_s (1%) [s]');
title('Tiempo de asentamiento al 1% vs coeficiente de amortiguamiento');
legend('t_s(c)', sprintf('c_c = %.2f N s/m', cc), 'Punto critico', 't_s minimo');
ylim([0 20]);

% Respuestas tipicas: sub, critico y sobreamortiguado
figure; hold on; grid on;
for c = [0.3*cc, cc, 1.8*cc]
    f = @(t,y) [y(2); -(c*y(2) + k*y(1))/m];
    [tt,y] = ode45(f, [0 5], [x0; v0], opts);
    plot(tt, y(:,1), 'LineWidth', 1.2);
end
plot([0 5], [tol tol], 'k:', [0 5], -[tol tol], 'k:');
xlabel('t [s]'); ylabel('x [m]');
legend('\zeta = 0.3', '\zeta = 1 (critico)', '\zeta = 1.8', 'banda 1%');
title('Respuesta libre para tres regimenes');
