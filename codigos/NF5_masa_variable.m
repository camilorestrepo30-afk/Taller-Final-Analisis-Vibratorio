%% NF5 - Oscilador con masa variable (recipiente que pierde arena)
% d/dt( m(t)*x' ) + c*x' + k*x = 0 ,  m(t) = M0 + mdot*t
% => m x'' + (mdot + c) x' + k x = 0
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
M0   = 10;     % masa inicial [kg]
mdot = -0.2;   % tasa de cambio de masa [kg/s]
mf   = 1;      % masa final [kg]
k    = 80;     % rigidez [N/m]
c    = 1.5;    % amortiguamiento [N s/m]
x0   = 0.1;    % desplazamiento inicial [m]
v0   = 0;
dt   = 1e-3;

tf = (mf - M0)/mdot;                % tiempo hasta que m(t) = 1 kg
fprintf('Tiempo de vaciado hasta m = %.1f kg: tf = %.2f s\n', mf, tf);

%% ---------- INTEGRACION ----------
mt = @(t) M0 + mdot*t;
f  = @(t,y) [y(2); -((mdot + c)*y(2) + k*y(1))/mt(t)];
t  = 0:dt:tf;
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~,y] = ode45(f, t, [x0; v0], opts);
x = y(:,1)';

%% ---------- FRECUENCIA INSTANTANEA ----------
% Numerica: cruces por cero ascendentes -> periodo de cada ciclo
iz = find(x(1:end-1) < 0 & x(2:end) >= 0);
tz = t(iz) - x(iz).*dt./(x(iz+1) - x(iz));     % interpolacion del cruce
Tn = diff(tz);
w_num = 2*pi./Tn;            t_num = tz(1:end-1) + Tn/2;
% Teorica (cuasi-estatica)
w_teo = sqrt(k./mt(t));

%% ---------- GRAFICAS ----------
figure;
subplot(3,1,1); plot(t, mt(t), 'k', 'LineWidth', 1.3); grid on;
ylabel('m(t) [kg]'); title('Masa del recipiente');
subplot(3,1,2); plot(t, x, 'b'); grid on;
ylabel('x [m]'); title('Desplazamiento');
subplot(3,1,3); plot(t_num, w_num, 'ro', t, w_teo, 'b-'); grid on;
xlabel('t [s]'); ylabel('\omega [rad/s]');
legend('Numerica (cruces por cero)', '\omega = sqrt(k/m(t))', 'Location', 'northwest');
title('Frecuencia instantanea');

fprintf('w inicial = %.3f rad/s , w final = %.3f rad/s\n', w_teo(1), w_teo(end));
