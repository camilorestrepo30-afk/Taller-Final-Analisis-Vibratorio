%% F1 - Barrido lineal de frecuencia (chirp) a traves de la resonancia
% m x'' + c x' + k x = F0 sen(w0 t + 1/2 beta t^2)  ->  w(t) = w0 + beta t
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
m  = 1;   k = 100;   c = 0.5;    % sistema
F0 = 1;                          % amplitud de la fuerza [N]
w0 = 2;                          % frecuencia inicial del barrido [rad/s]
w1 = 20;                         % frecuencia final del barrido [rad/s]
betas = [0.05 0.5 2 5];          % tasas de barrido [rad/s^2] (lenta -> rapida)
dt = 1e-3;

wn   = sqrt(k/m);  zeta = c/(2*sqrt(k*m));
wr   = wn*sqrt(1 - 2*zeta^2);                          % pico estacionario
Xest = (F0/k)/(2*zeta*sqrt(1 - zeta^2));               % amplitud max. estacionaria
fprintf('wn = %.3f rad/s , zeta = %.4f , w_pico = %.3f rad/s , Xmax est = %.4f m\n', ...
        wn, zeta, wr, Xest);

opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.01);
w_pico = zeros(size(betas));  X_pico = zeros(size(betas));
figure; hold on; grid on;
colores = lines(length(betas));

for i = 1:length(betas)
    b  = betas(i);
    tf = (w1 - w0)/b;                          % duracion del barrido
    t  = 0:dt:tf;
    F  = @(t) F0*sin(w0*t + 0.5*b*t.^2);
    f  = @(t,y) [y(2); (F(t) - c*y(2) - k*y(1))/m];
    [~,y] = ode45(f, t, [0; 0], opts);
    x  = y(:,1)';
    w_inst = w0 + b*t;                         % frecuencia instantanea

    [X_pico(i), ip] = max(abs(x));
    w_pico(i) = w_inst(ip);
    plot(w_inst, x, 'Color', colores(i,:));
    fprintf('beta = %5.2f rad/s^2 -> |x|max = %.4f m en w = %.3f rad/s (%.1f%% del estacionario)\n', ...
            b, X_pico(i), w_pico(i), 100*X_pico(i)/Xest);
end

% Respuesta estacionaria teorica (envolvente)
wv = linspace(w0, w1, 2000);
H  = (F0/k)./sqrt((1 - (wv/wn).^2).^2 + (2*zeta*wv/wn).^2);
plot(wv, H, 'k--', 'LineWidth', 1.5); plot(wv, -H, 'k--', 'LineWidth', 1.5);
xlabel('\omega(t) = \omega_0 + \beta t [rad/s]'); ylabel('x [m]');
legend([arrayfun(@(b) sprintf('\\beta = %.2f', b), betas, 'UniformOutput', false), {'|H| estacionaria'}]);
title('Respuesta al barrido de frecuencia');
xlim([5 15]);

figure;
subplot(2,1,1); semilogx(betas, w_pico, 'bo-'); hold on; grid on;
semilogx(betas, wr*ones(size(betas)), 'k--');
ylabel('\omega del pico [rad/s]'); title('Desplazamiento del pico con la tasa de barrido');
subplot(2,1,2); semilogx(betas, X_pico, 'ro-'); hold on; grid on;
semilogx(betas, Xest*ones(size(betas)), 'k--');
xlabel('\beta [rad/s^2]'); ylabel('|x|_{max} [m]');
