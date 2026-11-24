%% NF4 - Duffing libre: la frecuencia depende de la amplitud
% m*x'' + c*x' + k1*x + k3*x^3 = 0
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
m  = 1;      % masa [kg]
c  = 0.1;    % amortiguamiento [N s/m]
k1 = 10;     % rigidez lineal [N/m]
k3 = 50;     % rigidez cubica [N/m^3]
x0 = 1.5;    % desplazamiento inicial [m]  (grande para notar la no linealidad)
v0 = 0;
tf = 80;     % tiempo final [s]
dt = 1e-3;

%% ---------- INTEGRACION ----------
t    = 0:dt:tf;
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
fNL  = @(t,y) [y(2); -(c*y(2) + k1*y(1) + k3*y(1)^3)/m];   % no lineal
fL   = @(t,y) [y(2); -(c*y(2) + k1*y(1))/m];               % lineal (k3 = 0)
[~,yN] = ode45(fNL, t, [x0; v0], opts);   x  = yN(:,1)';
[~,yL] = ode45(fL,  t, [x0; v0], opts);   xL = yL(:,1)';

%% ---------- DETECTOR DE PICOS LOCALES (sin toolbox) ----------
p  = find(x(2:end-1)  > x(1:end-2)  & x(2:end-1)  >= x(3:end)  & x(2:end-1)  > 0) + 1;
pL = find(xL(2:end-1) > xL(1:end-2) & xL(2:end-1) >= xL(3:end) & xL(2:end-1) > 0) + 1;

% Periodo entre picos consecutivos -> frecuencia instantanea
T_NL = diff(t(p));        w_NL = 2*pi./T_NL;     A_NL = x(p(1:end-1));
T_L  = diff(t(pL));       w_L  = 2*pi./T_L;      A_L  = xL(pL(1:end-1));

% Aproximacion teorica (balance armonico de 1er orden)
A_teo = linspace(0, x0, 200);
w_teo = sqrt((k1 + 0.75*k3*A_teo.^2)/m);

%% ---------- TABLA ----------
fprintf('  Pico   Amplitud[m]   w_NL[rad/s]   w_lineal[rad/s]\n');
for i = 1:min(10,length(w_NL))
    fprintf('  %3d     %8.4f      %8.4f        %8.4f\n', i, A_NL(i), w_NL(i), w_L(min(i,end)));
end

%% ---------- GRAFICAS ----------
figure;
subplot(2,1,1);
plot(t, x, 'b'); hold on; plot(t(p), x(p), 'ro'); grid on;
xlabel('t [s]'); ylabel('x [m]'); title('Duffing libre y picos locales detectados');

subplot(2,1,2);
plot(A_NL, w_NL, 'bo', A_L, w_L, 'rs'); hold on; grid on;
plot(A_teo, w_teo, 'k-');
xlabel('Amplitud del pico [m]'); ylabel('\omega [rad/s]');
legend('No lineal (numerico)', 'Lineal (numerico)', '\omega \approx sqrt((k_1+3/4 k_3 A^2)/m)', ...
       'Location', 'northwest');
title('Frecuencia de oscilacion vs amplitud instantanea');
