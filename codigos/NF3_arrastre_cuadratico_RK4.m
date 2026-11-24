%% NF3 - Oscilador con arrastre cuadratico resuelto con RK4 propio
% m*x'' + cd*x'*|x'| + k*x = 0   vs   m*x'' + c_eq*x' + k*x = 0
clear; clc; close all;

%% ---------- GRAFICAS CON FONDO BLANCO ----------
set(groot, 'defaultFigureColor', 'w');          % fondo de la ventana
set(groot, 'defaultAxesColor', 'w');            % fondo de los ejes
set(groot, 'defaultAxesXColor', 'k', 'defaultAxesYColor', 'k', 'defaultAxesZColor', 'k');
set(groot, 'defaultTextColor', 'k', 'defaultAxesGridColor', [0.15 0.15 0.15]);
try, set(groot, 'defaultLegendColor', 'w', 'defaultLegendTextColor', 'k'); end
set(groot, 'defaultAxesFontSize', 11);
try, set(groot, 'defaultFigureTheme', 'light'); end   % MATLAB R2025a+ (modo oscuro)

%% ---------- PARAMETROS ----------
m  = 0.5;        % masa [kg]
k  = 20;         % rigidez [N/m]
cd = 0.8;        % coeficiente de arrastre cuadratico [kg/m]
x0 = 0.5;        % desplazamiento inicial [m]
v0 = 0;          % velocidad inicial [m/s]
tf = 20;         % tiempo final [s]
h  = 1e-3;       % paso del RK4 [s]

wn = sqrt(k/m);
% Amortiguador lineal equivalente (misma energia disipada en el 1er ciclo)
ceq = 8*cd*wn*x0/(3*pi);
fprintf('wn = %.4f rad/s ,  c_eq = %.4f N s/m\n', wn, ceq);

%% ---------- ECUACIONES (forma de estado y = [x; v]) ----------
f_cuad = @(t,y) [y(2); -(cd*y(2)*abs(y(2)) + k*y(1))/m];
f_lin  = @(t,y) [y(2); -(ceq*y(2) + k*y(1))/m];

%% ---------- RK4 DESDE CERO ----------
t  = 0:h:tf;  N = length(t);
yq = zeros(2,N);  yq(:,1) = [x0; v0];    % caso cuadratico
yl = zeros(2,N);  yl(:,1) = [x0; v0];    % caso lineal
for n = 1:N-1
    % --- cuadratico ---
    k1 = f_cuad(t(n),       yq(:,n));
    k2 = f_cuad(t(n)+h/2,   yq(:,n) + h/2*k1);
    k3 = f_cuad(t(n)+h/2,   yq(:,n) + h/2*k2);
    k4 = f_cuad(t(n)+h,     yq(:,n) + h*k3);
    yq(:,n+1) = yq(:,n) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    % --- lineal ---
    k1 = f_lin(t(n),        yl(:,n));
    k2 = f_lin(t(n)+h/2,    yl(:,n) + h/2*k1);
    k3 = f_lin(t(n)+h/2,    yl(:,n) + h/2*k2);
    k4 = f_lin(t(n)+h,      yl(:,n) + h*k3);
    yl(:,n+1) = yl(:,n) + h/6*(k1 + 2*k2 + 2*k3 + k4);
end
xq = yq(1,:);  xl = yl(1,:);

%% ---------- ENVOLVENTES ----------
% Picos numericos (maximos locales positivos)
pq = find(xq(2:end-1) > xq(1:end-2) & xq(2:end-1) > xq(3:end) & xq(2:end-1) > 0) + 1;
pl = find(xl(2:end-1) > xl(1:end-2) & xl(2:end-1) > xl(3:end) & xl(2:end-1) > 0) + 1;
% Envolventes teoricas
beta = 4*cd*wn/(3*pi*m);
env_q = x0 ./ (1 + beta*x0*t);         % decaimiento hiperbolico (cuadratico)
env_l = x0 * exp(-ceq*t/(2*m));        % decaimiento exponencial (lineal)

figure;
subplot(2,1,1);
plot(t, xq, 'b', t, xl, 'r'); grid on;
xlabel('t [s]'); ylabel('x [m]');
legend('Arrastre cuadratico', 'Amortiguamiento lineal equivalente');
title('Respuesta libre (RK4 propio)');

subplot(2,1,2);
semilogy(t(pq), xq(pq), 'bo', t(pl), xl(pl), 'rs'); hold on; grid on;
semilogy(t, env_q, 'b-', t, env_l, 'r-');
xlabel('t [s]'); ylabel('Amplitud de pico [m]');
legend('Picos cuadratico', 'Picos lineal', 'X_0/(1+\beta X_0 t)', 'X_0 e^{-c t/2m}');
title('Envolventes de decaimiento (escala logaritmica)');

fprintf('Ultimo pico (t = %.0f s): cuadratico = %.3e m , lineal = %.3e m\n', ...
        tf, xq(pq(end)), xl(pl(end)));
