%% F3 - Diagrama de bifurcacion del oscilador de Duffing (doble pozo)
% x'' + delta x' + alpha x + beta x^3 = gamma cos(w t)
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
delta = 0.3;  alpha = -1;  beta = 1;  w = 1.2;
gam   = linspace(0.20, 0.65, 1000);  % barrido del parametro de fuerza
n_trans = 100;                       % ciclos transitorios a descartar
n_reg   = 60;                        % ciclos registrados (puntos de Poincare)
ppc     = 200;                       % pasos de RK4 por periodo
x0 = 1;  v0 = 0;                     % condicion inicial (igual para todos los gamma)

T = 2*pi/w;   h = T/ppc;

%% ---------- RK4 VECTORIZADO (todos los gamma a la vez) ----------
% Cada columna corresponde a un valor de gamma -> mucho mas rapido que un bucle
fx = @(t,x,v) v;
fv = @(t,x,v) -delta*v - alpha*x - beta*x.^3 + gam*cos(w*t);

x = x0*ones(size(gam));  v = v0*ones(size(gam));  t = 0;
P = zeros(n_reg, length(gam));       % secciones de Poincare

for ciclo = 1:(n_trans + n_reg)
    for s = 1:ppc
        k1x = fx(t,x,v);                    k1v = fv(t,x,v);
        k2x = fx(t+h/2,x+h/2*k1x,v+h/2*k1v); k2v = fv(t+h/2,x+h/2*k1x,v+h/2*k1v);
        k3x = fx(t+h/2,x+h/2*k2x,v+h/2*k2v); k3v = fv(t+h/2,x+h/2*k2x,v+h/2*k2v);
        k4x = fx(t+h,x+h*k3x,v+h*k3v);       k4v = fv(t+h,x+h*k3x,v+h*k3v);
        x = x + h/6*(k1x + 2*k2x + 2*k3x + k4x);
        v = v + h/6*(k1v + 2*k2v + 2*k3v + k4v);
        t = t + h;
    end
    if ciclo > n_trans
        P(ciclo - n_trans, :) = x;     % muestreo estroboscopico cada T = 2*pi/w
    end
end

%% ---------- DIAGRAMA DE BIFURCACION ----------
figure;
G = repmat(gam, n_reg, 1);
plot(G(:), P(:), 'k.', 'MarkerSize', 2); grid on;
xlabel('\gamma'); ylabel('x(nT)');
title('Diagrama de bifurcacion (seccion de Poincare) - Duffing');

% Numero de puntos distintos por gamma: 1 = periodo-1, 2 = periodo-2, muchos = caos
nper = zeros(size(gam));
for j = 1:length(gam)
    nper(j) = length(unique(round(P(:,j)*1e3)/1e3));
end
for gsel = [0.20 0.28 0.29 0.37 0.50 0.65]
    [~,j] = min(abs(gam - gsel));
    fprintf('gamma = %.3f -> %2d valores distintos en Poincare\n', gam(j), nper(j));
end

%% ---------- SECCION DE POINCARE (x, v) PARA UN GAMMA CAOTICO ----------
gc  = 0.50;                          % gamma a estudiar (modificable)
Nc  = 2000;                          % numero de cruces de Poincare
fgc = @(t,y) [y(2); -delta*y(2) - alpha*y(1) - beta*y(1)^3 + gc*cos(w*t)];
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
[~,y] = ode45(fgc, (0:Nc+n_trans)*T, [x0; v0], opts);   % salida solo cada T
y = y(n_trans+1:end, :);
figure; plot(y(:,1), y(:,2), 'k.', 'MarkerSize', 3); grid on;
xlabel('x(nT)'); ylabel('v(nT)');
title(sprintf('Seccion de Poincare, \\gamma = %.2f (atractor extrano)', gc));
