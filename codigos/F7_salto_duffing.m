%% F7 - Fenomeno de salto en Duffing endurecedor (barrido ascendente y descendente)
% x'' + 0.1 x' + x + 0.5 x^3 = 1.5 cos(w t)
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
c = 0.1;  k1 = 1;  k3 = 0.5;  F0 = 1.5;
w_min = 0.5;  w_max = 5.0;  Nw = 120;    % numero de frecuencias
n_ciclos = 60;                           % ciclos simulados por frecuencia
n_med    = 10;                           % ultimos ciclos para medir amplitud
ppc      = 100;                          % muestras por ciclo

w_sube = linspace(w_min, w_max, Nw);
w_baja = fliplr(w_sube);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

A_sube = zeros(1,Nw);  A_baja = zeros(1,Nw);

%% ---------- BARRIDO ASCENDENTE ----------
y0 = [0; 0];
for i = 1:Nw
    w = w_sube(i);  T = 2*pi/w;
    f = @(t,y) [y(2); F0*cos(w*t) - c*y(2) - k1*y(1) - k3*y(1)^3];
    t = linspace(0, n_ciclos*T, n_ciclos*ppc + 1);
    [~,y] = ode45(f, t, y0, opts);
    xs = y(end - n_med*ppc:end, 1);
    A_sube(i) = (max(xs) - min(xs))/2;
    y0 = y(end,:)';          % continuacion: estado final -> condicion inicial
end

%% ---------- BARRIDO DESCENDENTE ----------
% (y0 continua desde el final del barrido ascendente)
for i = 1:Nw
    w = w_baja(i);  T = 2*pi/w;
    f = @(t,y) [y(2); F0*cos(w*t) - c*y(2) - k1*y(1) - k3*y(1)^3];
    t = linspace(0, n_ciclos*T, n_ciclos*ppc + 1);
    [~,y] = ode45(f, t, y0, opts);
    xs = y(end - n_med*ppc:end, 1);
    A_baja(i) = (max(xs) - min(xs))/2;
    y0 = y(end,:)';
end

%% ---------- SALTOS ----------
[~,is] = min(diff(A_sube));   % mayor caida al subir
[~,ib] = max(diff(A_baja));
fprintf('Salto hacia abajo (subiendo w): entre %.3f y %.3f rad/s\n', w_sube(is), w_sube(is+1));
fprintf('Salto hacia arriba (bajando w): entre %.3f y %.3f rad/s\n', w_baja(ib), w_baja(ib+1));
fprintf('Zona biestable aprox.: %.2f < w < %.2f rad/s\n', w_baja(ib+1), w_sube(is));

%% ---------- CURVA TEORICA (balance armonico) ----------
% [(k1 - w^2 + 3/4 k3 A^2)^2 + (c w)^2] A^2 = F0^2
Av = linspace(0.01, 6, 3000);  wt = []; At = [];
for A = Av
    % resolver en w^2: (a - W)^2 + c^2 W = (F0/A)^2 , a = k1 + 0.75 k3 A^2
    a = k1 + 0.75*k3*A^2;
    r = roots([1, c^2 - 2*a, a^2 - (F0/A)^2]);
    r = r(abs(imag(r)) < 1e-9 & real(r) > 0);
    wt = [wt; sqrt(real(r))];  At = [At; A*ones(length(r),1)];
end

figure;
plot(w_sube, A_sube, 'b>-', w_baja, A_baja, 'r<-'); hold on; grid on;
plot(wt, At, 'k.', 'MarkerSize', 3);
xlabel('\omega [rad/s]'); ylabel('Amplitud estacionaria [m]');
legend('\omega creciente', '\omega decreciente', 'Balance armonico');
title('Curvas de respuesta en frecuencia: salto e histeresis');
xlim([w_min w_max]);
