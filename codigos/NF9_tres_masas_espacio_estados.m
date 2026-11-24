%% NF9 - Cadena de 3 masas con amortiguamiento NO proporcional (espacio de estados 6x6)
% pared --k1,c1-- m1 --k2,c2-- m2 --k3,c3-- m3 --k4,c4-- pared
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
m  = 1;                          % masa de cada bloque [kg]
kv = [100 150 200 120];          % rigideces k1..k4 [N/m]
cv = [0.5 3.0 0.1 1.0];          % amortiguadores c1..c4 [N s/m] (distintos)
x0 = [0.1; 0; 0];                % desplazamiento inicial SOLO en la masa 1 [m]
v0 = [0; 0; 0];
tf = 30;  dt = 1e-3;

%% ---------- MATRICES M, K, C ----------
M = m*eye(3);
K = [kv(1)+kv(2), -kv(2),       0;
     -kv(2),      kv(2)+kv(3), -kv(3);
     0,          -kv(3),        kv(3)+kv(4)];
C = [cv(1)+cv(2), -cv(2),       0;
     -cv(2),      cv(2)+cv(3), -cv(3);
     0,          -cv(3),        cv(3)+cv(4)];

% Verificacion: C NO es proporcional (condicion de Caughey)
cau = norm(C/M*K - K/M*C);
fprintf('||C M^-1 K - K M^-1 C|| = %.4f  (distinto de 0 -> no proporcional)\n', cau);

%% ---------- ESPACIO DE ESTADOS z = [x; v] ----------
A = [zeros(3), eye(3);
     -M\K,     -M\C];
disp('Matriz A (6x6):'); disp(A);

%% ---------- AUTOVALORES COMPLEJOS ----------
lam = eig(A);
lam = lam(imag(lam) > 0);              % un valor de cada par conjugado
[~,ord] = sort(abs(lam)); lam = lam(ord);
wn   = abs(lam);                       % frecuencia natural
zeta = -real(lam)./abs(lam);           % razon de amortiguamiento
wd   = imag(lam);                      % frecuencia amortiguada
fprintf('\n Modo   Re(lam)     Im(lam)=wd    wn[rad/s]   zeta\n');
for i = 1:3
    fprintf('  %d   %9.4f   %9.4f    %9.4f   %7.4f\n', i, real(lam(i)), wd(i), wn(i), zeta(i));
end

%% ---------- RESPUESTA LIBRE ----------
t = 0:dt:tf;
z0 = [x0; v0];
opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
[~,z] = ode45(@(t,z) A*z, t, z0, opts);

% Comprobacion con la solucion exacta z(t) = expm(A t) z0
z_ex = expm(A*t(end))*z0;
fprintf('\nError ode45 vs expm en t = %.0f s: %.2e\n', tf, norm(z(end,:)' - z_ex));

figure;
for i = 1:3
    subplot(3,1,i); plot(t, z(:,i)); grid on;
    ylabel(sprintf('x_%d [m]', i));
end
xlabel('t [s]'); subplot(3,1,1); title('Respuesta libre con x_1(0) = 0.1 m');

figure;
plot(real(eig(A)), imag(eig(A)), 'rx', 'MarkerSize', 10, 'LineWidth', 2); grid on;
xlabel('Re(\lambda)'); ylabel('Im(\lambda)'); title('Autovalores en el plano complejo');
